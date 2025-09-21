#!/usr/bin/env python3
"""Run MSFragger + Philosopher on a local dataset and emit protein matrices.

This utility converts (if necessary) mzML.gz inputs, executes MSFragger with the
provided parameter file, runs the Philosopher workflow (peptideprophet ->
proteinprophet -> filter -> report), and generates convenient downstream
artifacts:

* `protein.tsv` / `peptide.tsv` / `psm.tsv` / `ion.tsv` from Philosopher
* `protein_abundance.tsv` – log2-intensity matrix indexed by Gene/Protein
* `metadata.tsv` – sample annotations parsed from filenames
* `uniprot_to_hgnc.tsv` – UniProt -> gene mapping used for collapsing IDs

Example
-------
```bash
micromamba run -n proteomics python scripts/py/proteomics_pipeline.py \
  --pxd PXD011796 \
  --mzml-dir data/interim/proteomics/PXD011796/mzML \
  --workspace data/interim/proteomics/PXD011796/fragger_closed \
  --params resources/msfragger/default_closed.params \
  --fasta data/interim/proteomics/reference/UP000005640_9606_td.fasta \
  --output-dir data/processed/proteomics/PXD011796
```
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


DEFAULT_MSFRAGGER_JAR = Path("env/.mamba/envs/proteomics/share/msfragger-4.1-0/MSFragger.jar")
DEFAULT_PHILOSOPHER_BIN = Path("env/bin/philosopher")
UNIPROT_HEADERS = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/118.0.0.0 Safari/537.36 nets-fuction/1.0 (+https://example.org/contact)",
    "Accept": "application/json",
}
UNIPROT_MAX_BATCH = 100
UNIPROT_CHUNK_DELAY = 0.3
UNIPROT_MAX_RETRIES = 4
UNIPROT_RETRY_BACKOFF = 0.7

logger = logging.getLogger("proteomics_pipeline")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def run_command(cmd: List[str], cwd: Optional[Path] = None) -> None:
    logger.debug("Running command: %s", " ".join(cmd))
    proc = subprocess.run(cmd, cwd=cwd, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}")


def ensure_mzml_inputs(mzml_dir: Path, workspace: Path, limit: Optional[int] = None) -> List[Path]:
    """Collect mzML files, decompressing gz archives into workspace if needed."""

    if not mzml_dir.exists():
        raise FileNotFoundError(f"mzML directory not found: {mzml_dir}")

    mzml_inputs: List[Path] = []
    target_dir = workspace / "mzml"
    target_dir.mkdir(parents=True, exist_ok=True)

    for path in sorted(mzml_dir.glob("**/*")):
        if path.is_dir():
            continue
        suffix = path.suffix.lower()
        if suffix == ".mzml":
            mzml_inputs.append(path)
        elif suffix == ".gz" and path.name.lower().endswith(".mzml.gz"):
            dest = target_dir / path.stem  # removes .gz
            if not dest.exists() or dest.stat().st_mtime < path.stat().st_mtime:
                logger.info("Decompressing %s -> %s", path.name, dest)
                with gzip.open(path, "rb") as src, open(dest, "wb") as out:
                    shutil.copyfileobj(src, out)
            else:
                logger.debug("Using cached decompressed file %s", dest)
            mzml_inputs.append(dest)
        if limit is not None and len(mzml_inputs) >= limit:
            break
    if limit is not None:
        mzml_inputs = mzml_inputs[:limit]

    if not mzml_inputs:
        raise FileNotFoundError(f"No mzML or mzML.gz files found under {mzml_dir}")
    return mzml_inputs


def run_msfragger(msfragger_jar: Path, params: Path, mzml_files: Iterable[Path], memory_gb: int, workspace: Path) -> None:
    cmd = [
        "java",
        f"-Xmx{memory_gb}G",
        "-jar",
        str(msfragger_jar),
        str(params),
    ] + [str(p) for p in mzml_files]
    run_command(cmd, cwd=workspace)


def run_philosopher(philosopher_bin: Path, fasta: Path, workspace: Path) -> None:
    run_command([str(philosopher_bin), "workspace", "--clean"], cwd=workspace)
    run_command([str(philosopher_bin), "workspace", "--init"], cwd=workspace)
    run_command([str(philosopher_bin), "database", "--custom", str(fasta)], cwd=workspace)
    run_command([str(philosopher_bin), "database", "--annotate", str(fasta)], cwd=workspace)

    pepxmls = sorted(workspace.rglob("*.pepXML"))
    if not pepxmls:
        raise FileNotFoundError("No pepXML files produced by MSFragger")

    local_pepxmls: List[Path] = []
    for src in pepxmls:
        dest = workspace / src.name
        if src.parent != workspace:
            shutil.copy(src, dest)
        else:
            dest = src
        local_pepxmls.append(dest)

    run_command(
        [
            str(philosopher_bin),
            "peptideprophet",
            "--nonparam",
            "--ppm",
            "--expectscore",
            "--decoy",
            "rev_",
        ]
        + [str(p) for p in local_pepxmls],
        cwd=workspace,
    )

    interact = sorted(workspace.glob("interact-*.pep.xml"))
    run_command([str(philosopher_bin), "proteinprophet", "--maxppmdiff", "200", "--output", "combined"] + [str(f) for f in interact], cwd=workspace)
    filter_cmd = [
        str(philosopher_bin),
        "filter",
        "--sequential",
        "--razor",
        "--mapmods",
        "--tag",
        "rev_",
    ]
    for f in interact:
        filter_cmd.extend(["--pepxml", str(f.name)])
    filter_cmd.extend(["--protxml", "combined.prot.xml"])
    run_command(filter_cmd, cwd=workspace)
    run_command([str(philosopher_bin), "report", "--msstats"], cwd=workspace)


def load_protein_report(workspace: Path) -> pd.DataFrame:
    candidates = sorted(workspace.glob("*protein.tsv"))
    if not candidates:
        raise FileNotFoundError("Philosopher protein report not found")
    logger.info("Using protein report: %s", candidates[0].name)
    return pd.read_csv(candidates[0], sep="\t")


def _create_uniprot_session() -> requests.Session:
    """Instantiate a requests session tuned for UniProt API access."""

    retry_kwargs = dict(
        total=UNIPROT_MAX_RETRIES,
        read=UNIPROT_MAX_RETRIES,
        connect=UNIPROT_MAX_RETRIES,
        status=UNIPROT_MAX_RETRIES,
        backoff_factor=UNIPROT_RETRY_BACKOFF,
        status_forcelist=(429, 500, 502, 503, 504),
    )
    try:
        retry = Retry(allowed_methods=("GET", "POST"), **retry_kwargs)
    except TypeError:  # urllib3 < 1.26 uses method_whitelist
        retry = Retry(method_whitelist=("GET", "POST"), **retry_kwargs)
    adapter = HTTPAdapter(max_retries=retry)
    session = requests.Session()
    session.headers.update(UNIPROT_HEADERS)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def map_uniprot_to_gene(uniprot_ids: List[str]) -> pd.DataFrame:
    ids = [uid for uid in set(uniprot_ids) if isinstance(uid, str) and uid]
    if not ids:
        return pd.DataFrame(columns=["UniProt", "Gene"])

    logger.info("Submitting UniProt mapping for %d unique accessions", len(ids))
    session = _create_uniprot_session()
    rows: List[Dict[str, str]] = []

    for start in range(0, len(ids), UNIPROT_MAX_BATCH):
        chunk = ids[start : start + UNIPROT_MAX_BATCH]
        try:
            rows.extend(_map_uniprot_chunk(chunk, session))
        except Exception as exc:
            logger.warning("UniProt mapping chunk failed (%s entries): %s", len(chunk), exc)
        time.sleep(UNIPROT_CHUNK_DELAY)

    if not rows:
        logger.warning("UniProt mapping returned no records; falling back to accession IDs")
        return pd.DataFrame(columns=["UniProt", "Gene"])

    return pd.DataFrame(rows).drop_duplicates()


def _map_uniprot_chunk(chunk: List[str], session: requests.Session) -> List[Dict[str, str]]:
    attempt = 0
    while True:
        attempt += 1
        try:
            response = session.post(
                "https://rest.uniprot.org/idmapping/run",
                data={"from": "UniProtKB_AC-ID", "to": "Gene_Name", "ids": ",".join(chunk)},
                timeout=30,
            )
            if response.status_code == 403:
                raise requests.HTTPError("UniProt request rejected with 403", response=response)
            response.raise_for_status()
            payload = response.json()
            job_id = payload.get("jobId")
            if not job_id:
                raise RuntimeError("UniProt job submission failed: missing jobId")

            status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
            result_url = f"https://rest.uniprot.org/idmapping/stream/{job_id}"
            while True:
                status = session.get(status_url, timeout=30)
                if status.status_code == 403:
                    raise requests.HTTPError("UniProt status poll rejected with 403", response=status)
                status.raise_for_status()
                state = status.json().get("jobStatus")
                if state == "FINISHED":
                    break
                if state == "FAILED":
                    raise RuntimeError("UniProt mapping failed")
                time.sleep(1)

            results = session.get(result_url, timeout=30)
            if results.status_code == 403:
                raise requests.HTTPError("UniProt result fetch rejected with 403", response=results)
            results.raise_for_status()
            data = results.json().get("results", [])
            chunk_rows: List[Dict[str, str]] = []
            for entry in data:
                gene = entry.get("to")
                if not gene:
                    continue
                chunk_rows.append({"UniProt": entry["from"], "Gene": gene.split(";")[0]})
            return chunk_rows
        except requests.RequestException as exc:
            if attempt >= UNIPROT_MAX_RETRIES:
                raise
            delay = UNIPROT_RETRY_BACKOFF * (2 ** (attempt - 1))
            logger.debug(
                "Retrying UniProt chunk (attempt %s/%s, size %s) after %.1fs due to %s",
                attempt,
                UNIPROT_MAX_RETRIES,
                len(chunk),
                delay,
                exc,
            )
            time.sleep(delay)


def build_protein_matrix(df: pd.DataFrame, output_dir: Path, skip_uniprot: bool = False) -> pd.DataFrame:
    intensity_cols = [c for c in df.columns if "Intensity" in c]
    if not intensity_cols:
        raise ValueError("No intensity columns detected in protein report")

    df = df.copy()
    df["PrimaryAcc"] = df["Protein"].astype(str).str.split(";").str[0]
    mapping = pd.DataFrame()
    if not skip_uniprot:
        try:
            mapping = map_uniprot_to_gene(df["PrimaryAcc"].tolist())
        except Exception as exc:
            logger.warning("UniProt mapping failed: %s", exc)
            mapping = pd.DataFrame()

    if not mapping.empty:
        df = df.merge(mapping, how="left", left_on="PrimaryAcc", right_on="UniProt")
        df["Gene"] = df["Gene"].fillna(df["PrimaryAcc"])
    else:
        logger.warning("UniProt mapping empty; falling back to accession IDs")
        df["Gene"] = df["PrimaryAcc"]

    df = df.drop_duplicates(subset=["Gene"]).set_index("Gene")
    matrix = np.log2(df[intensity_cols].replace(0, np.nan))

    matrix_path = output_dir / "protein_abundance.tsv"
    matrix.to_csv(matrix_path, sep="\t")
    logger.info("Protein matrix written to %s", matrix_path)

    map_path = output_dir / "uniprot_to_hgnc.tsv"
    if not mapping.empty:
        mapping.to_csv(map_path, sep="\t", index=False)
    else:
        df[["PrimaryAcc"]].rename(columns={"PrimaryAcc": "UniProt"}).assign(Gene=df.index).to_csv(map_path, sep="\t", index=False)

    return matrix


def derive_metadata(matrix: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    records = []
    for sample in matrix.columns:
        clean = sample.replace(".raw", "").replace(".mzML", "").replace(".pepXML", "")
        tokens = clean.split("_")
        disease = "Unknown"
        stimulus = "NA"
        timepoint = "NA"
        for token in tokens:
            upper = token.upper()
            if upper.startswith("SLE"):
                disease = "SLE"
            elif upper.startswith("RA"):
                disease = "RA"
            elif upper.startswith("NC"):
                disease = "Healthy"
            if upper in {"PMA", "A23", "A23187"}:
                stimulus = upper
            if upper.endswith("H") and upper[:-1].isdigit():
                timepoint = upper
        records.append({
            "sample_id": sample,
            "disease": disease,
            "stimulus": stimulus,
            "timepoint": timepoint,
        })
    meta = pd.DataFrame(records)
    meta_path = output_dir / "metadata.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)
    logger.info("Metadata written to %s", meta_path)
    return meta


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--pxd", required=True, help="Dataset identifier (used for logging only)")
    parser.add_argument("--mzml-dir", type=Path, required=True, help="Directory containing mzML or mzML.gz files")
    parser.add_argument("--workspace", type=Path, required=True, help="Working directory for MSFragger/Philosopher outputs")
    parser.add_argument("--params", type=Path, required=True, help="MSFragger parameter file")
    parser.add_argument("--fasta", type=Path, required=True, help="FASTA file (target-decoy) for MSFragger/Philosopher")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory to store processed matrices and reports")
    parser.add_argument("--msfragger-jar", type=Path, default=DEFAULT_MSFRAGGER_JAR, help="Path to MSFragger.jar")
    parser.add_argument("--philosopher-bin", type=Path, default=DEFAULT_PHILOSOPHER_BIN, help="Philosopher executable path")
    parser.add_argument("--memory-gb", type=int, default=32, help="Java heap size for MSFragger")
    parser.add_argument("--max-files", type=int, default=None, help="Process only the first N mzML files (for smoke tests)")
    parser.add_argument("--skip-uniprot", action="store_true", help="Skip UniProt mapping (use accession IDs)")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging level")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(levelname)s: %(message)s")

    if not args.msfragger_jar.exists():
        raise FileNotFoundError(f"MSFragger jar not found: {args.msfragger_jar}")
    if not args.philosopher_bin.exists():
        raise FileNotFoundError(f"Philosopher binary not found: {args.philosopher_bin}")
    if not args.params.exists():
        raise FileNotFoundError(f"Parameter file not found: {args.params}")
    if not args.fasta.exists():
        raise FileNotFoundError(f"FASTA file not found: {args.fasta}")

    workspace = args.workspace.resolve()
    output_dir = args.output_dir.resolve()
    mzml_dir = args.mzml_dir.resolve()

    workspace.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Processing %s", args.pxd)

    mzml_files = ensure_mzml_inputs(mzml_dir, workspace, args.max_files)
    logger.info("Found %d mzML inputs (max=%s)", len(mzml_files), args.max_files)

    msfragger_jar = args.msfragger_jar.resolve()
    params = args.params.resolve()
    fasta = args.fasta.resolve()
    philosopher_bin = args.philosopher_bin.resolve()

    run_msfragger(msfragger_jar, params, mzml_files, args.memory_gb, workspace)
    run_philosopher(philosopher_bin, fasta, workspace)

    protein_df = load_protein_report(workspace)
    matrix = build_protein_matrix(protein_df, output_dir, skip_uniprot=args.skip_uniprot)
    derive_metadata(matrix, output_dir)

    logger.info("Completed pipeline for %s", args.pxd)


if __name__ == "__main__":
    try:
        import requests  # noqa: F401
    except ImportError:
        sys.stderr.write("ERROR: Python 'requests' package required. Install in proteomics env.\n")
        sys.exit(1)
    main()
