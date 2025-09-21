#!/usr/bin/env python3
"""Metabolomics Workbench archive extraction and redox ratio computation.

This script mimics the original Colab notebook logic but targets the local
project layout. It extracts each metabolomics archive into
``data/interim/metabolomics/<study>/raw`` (or a user-specified directory),
tidies the tables, harmonises a handful of redox-related identifiers, and
emits a wide table of ratios per sample.

Example
-------
```bash
scripts/py/metabolomics_preprocessing.py \
  --raw-dir data/raw/metabolomics \
  --extract-root data/interim/metabolomics/raw_tables \
  --output data/processed/metabolomics/redox_ratios.tsv
```
"""

from __future__ import annotations

import argparse
import logging
import re
import tarfile
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd

from mw_parser import download_named_metabolites_via_api, to_redox_canonical


STUDY_SYNONYMS: Dict[str, set[str]] = {
    "GSH": {"gsh", "reduced glutathione", "glutathione (reduced)", "hmdb0000125"},
    "GSSG": {"gssg", "oxidized glutathione", "glutathione disulfide", "hmdb0003336"},
    "NADPH": {
        "nadph",
        "beta-nicotinamide adenine dinucleotide phosphate reduced",
        "hmdb0000221",
    },
    "NADP": {"nadp", "nadp+", "hmdb0000219"},
    "MetSO": {"methionine sulfoxide", "l-methionine sulfoxide", "hmdb0003337"},
    "Met": {"methionine", "l-methionine", "hmdb0000696"},
    "Lactate": {"lactate", "l-lactic acid", "hmdb0000190"},
    "Pyruvate": {"pyruvate", "pyruvic acid", "hmdb0000243"},
}

RATIOS: Dict[str, tuple[str, str]] = {
    "GSH_GSSG": ("GSH", "GSSG"),
    "NADPH_NADP": ("NADPH", "NADP"),
    "MetSO_Met": ("MetSO", "Met"),
    "Lactate_Pyruvate": ("Lactate", "Pyruvate"),
}

STUDY_ID_PATTERN = re.compile(r"(ST\d{6})", re.IGNORECASE)


def safe_extract(archive: Path, dest: Path, force: bool = False) -> None:
    """Extract *archive* into *dest* if missing or if ``force`` is True."""

    dest.mkdir(parents=True, exist_ok=True)
    if any(dest.iterdir()) and not force:
        logging.debug("skip extraction for %s (already populated)", archive)
        return

    suffix = archive.suffix.lower()
    if suffix == ".zip":
        with zipfile.ZipFile(archive, "r") as zf:
            zf.extractall(dest)
    elif suffix in {".tar", ".gz", ".bz2", ".tgz"} or archive.name.endswith(
        (".tar.gz", ".tar.bz2")
    ):
        with tarfile.open(archive, "r:*") as tf:
            tf.extractall(dest)
    else:
        raise ValueError(f"Unsupported archive format: {archive}")


def find_tables(root: Path, patterns: Sequence[str] | None = None) -> List[Path]:
    pats = patterns or ["*.csv", "*.tsv", "*.txt", "*.xlsx"]
    files: List[Path] = []
    for pattern in pats:
        files.extend(root.rglob(pattern))
    return files


def infer_study_id(dataset_name: str) -> Optional[str]:
    match = STUDY_ID_PATTERN.search(dataset_name)
    if match:
        return match.group(1).upper()
    return None


def load_table(path: Path) -> pd.DataFrame:
    try:
        suffix = path.suffix.lower()
        if suffix == ".csv":
            return pd.read_csv(path)
        if suffix in {".tsv", ".txt"}:
            return pd.read_csv(path, sep="\t")
        if suffix == ".xlsx":
            return pd.read_excel(path)
    except Exception as exc:  # pragma: no cover - defensive log
        logging.warning("failed to load %s: %s", path, exc)
    return pd.DataFrame()


def _normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(col).strip() for col in df.columns]
    return df


def _detect_hmdb_column(df: pd.DataFrame) -> str | None:
    for col in df.columns:
        if "hmdb" in col.lower():
            return col
    return None


def _detect_metabolite_column(df: pd.DataFrame) -> str | None:
    for col in df.columns:
        lower = col.lower()
        if any(keyword in lower for keyword in ("metabolite", "compound", "analyte")):
            return col
        if lower == "name":
            return col
    return None


def melt_table(df: pd.DataFrame) -> pd.DataFrame:
    df = _normalise_columns(df)
    id_col = _detect_hmdb_column(df) or _detect_metabolite_column(df)
    if id_col is None:
        raise ValueError("no HMDB/metabolite column detected")

    value_cols = [c for c in df.columns if c != id_col]
    tidy = df.melt(id_vars=id_col, value_vars=value_cols, var_name="sample", value_name="intensity")
    tidy = tidy.dropna(subset=["intensity"])
    tidy["feature_id"] = tidy[id_col].astype(str).str.strip()
    tidy["sample"] = tidy["sample"].astype(str).str.strip()
    tidy = tidy[tidy["feature_id"] != ""]
    tidy = tidy.replace({"intensity": {"": np.nan}}).dropna(subset=["intensity"])
    tidy["intensity"] = pd.to_numeric(tidy["intensity"], errors="coerce")
    tidy = tidy.dropna(subset=["intensity"])
    return tidy[["feature_id", "sample", "intensity"]]


def harmonise_feature_id(feature: str) -> str:
    lowered = feature.lower()
    for label, synonyms in STUDY_SYNONYMS.items():
        if lowered in synonyms:
            return label
    canonical = to_redox_canonical(feature)
    if canonical:
        return canonical
    return feature


def calculate_ratios(tidy: pd.DataFrame) -> pd.DataFrame:
    tidy = tidy.copy()
    tidy["harmonised_id"] = tidy["feature_id"].apply(harmonise_feature_id)

    ratio_frames: List[pd.DataFrame] = []
    for ratio_name, (numerator, denominator) in RATIOS.items():
        num = tidy[tidy["harmonised_id"].str.lower() == numerator.lower()]
        den = tidy[tidy["harmonised_id"].str.lower() == denominator.lower()]
        if num.empty or den.empty:
            continue
        merged = num.merge(den, on="sample", suffixes=("_num", "_den"))
        merged = merged[merged["intensity_den"] != 0]
        if merged.empty:
            continue
        merged["metric"] = ratio_name
        merged["value"] = merged["intensity_num"] / merged["intensity_den"]
        merged = merged.replace({"value": [np.inf, -np.inf]}, np.nan).dropna(subset=["value"])
        ratio_frames.append(merged[["sample", "metric", "value"]])

    if not ratio_frames:
        return pd.DataFrame()

    combined = pd.concat(ratio_frames, ignore_index=True)
    pivot = combined.pivot_table(index="sample", columns="metric", values="value", aggfunc="mean")
    pivot.index.name = "sample_id"
    return pivot


def calculate_ratios_from_matrix(matrix: pd.DataFrame) -> pd.DataFrame:
    """Compute redox ratios from a named-metabolite matrix (metabolites × samples)."""

    if matrix is None or matrix.empty:
        return pd.DataFrame()

    df = matrix.copy()
    df.columns = [str(c).strip() for c in df.columns]

    candidate_cols = [c for c in df.columns if re.search(r"metabolite|compound|name", c.lower())]
    metabolite_col = candidate_cols[0] if candidate_cols else df.columns[0]

    # Normalise metabolite identifiers
    df = df[df[metabolite_col].notna()]
    df[metabolite_col] = df[metabolite_col].astype(str).str.strip()

    # Convert numeric columns; drop columns with no numeric observations
    numeric_cols: Dict[str, pd.Series] = {}
    for col in df.columns:
        if col == metabolite_col:
            continue
        series = pd.to_numeric(df[col], errors="coerce")
        if series.notna().any():
            numeric_cols[col] = series

    if not numeric_cols:
        return pd.DataFrame()

    numeric_df = pd.DataFrame(numeric_cols)
    numeric_df = numeric_df.replace({np.inf: np.nan, -np.inf: np.nan})
    tidy_df = pd.concat([df[[metabolite_col]].reset_index(drop=True), numeric_df.reset_index(drop=True)], axis=1)
    tidy_df = tidy_df.dropna(how="all", subset=numeric_df.columns.tolist())

    long = tidy_df.melt(id_vars=metabolite_col, var_name="sample", value_name="intensity")
    long = long.dropna(subset=["intensity"])
    long["feature_id"] = long[metabolite_col].astype(str).str.strip()
    long["sample"] = long["sample"].astype(str).str.strip()

    tidy = long[["feature_id", "sample", "intensity"]]
    tidy["intensity"] = pd.to_numeric(tidy["intensity"], errors="coerce")
    tidy = tidy.dropna(subset=["intensity"])

    return calculate_ratios(tidy)


def process_archive(archive: Path, extract_root: Path, force_extract: bool = False) -> pd.DataFrame:
    dataset_name = archive.stem

    study_id = infer_study_id(dataset_name)
    if study_id:
        logging.info("attempting REST datatable fetch for %s", study_id)
        try:
            matrix = download_named_metabolites_via_api(study_id)
        except Exception as exc:  # pragma: no cover - REST fallback
            logging.warning("REST fetch failed for %s: %s", study_id, exc)
            matrix = None
        if matrix is not None and not matrix.empty:
            source = matrix.attrs.get("source")
            logging.info("received named metabolite table via REST (%s)", source or "unknown source")
            ratios = calculate_ratios_from_matrix(matrix)
            if not ratios.empty:
                ratios["dataset"] = dataset_name
                return ratios
            logging.warning("REST table for %s did not contain required redox metabolites", dataset_name)
        else:
            logging.info("no REST datatable available for %s", study_id)

    extract_dir = extract_root / dataset_name
    logging.info("extracting %s -> %s", archive.name, extract_dir)
    safe_extract(archive, extract_dir, force=force_extract)

    ratio_frames: List[pd.DataFrame] = []
    for table_path in find_tables(extract_dir):
        df = load_table(table_path)
        if df.empty:
            continue
        try:
            tidy = melt_table(df)
        except ValueError:
            continue
        ratios = calculate_ratios(tidy)
        if ratios.empty:
            continue
        ratios["dataset"] = dataset_name
        ratio_frames.append(ratios)

    if not ratio_frames:
        logging.warning("no ratios generated for %s", archive.name)
        return pd.DataFrame()

    merged = pd.concat(ratio_frames)
    merged = merged.groupby(level=0).mean()
    merged["dataset"] = dataset_name
    return merged


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=Path("data/raw/metabolomics"),
        help="Directory containing metabolomics archives (default: data/raw/metabolomics)",
    )
    parser.add_argument(
        "--extract-root",
        type=Path,
        default=Path("data/interim/metabolomics/raw_tables"),
        help="Root directory for extracted archives",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/metabolomics/redox_ratios.tsv"),
        help="Output TSV path for the redox ratio table",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="*.zip",
        help="Glob pattern for archives (default: *.zip)",
    )
    parser.add_argument(
        "--force-extract",
        action="store_true",
        help="Re-extract archives even if a non-empty directory already exists",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(levelname)s: %(message)s")

    raw_dir: Path = args.raw_dir
    extract_root: Path = args.extract_root
    output_path: Path = args.output

    extract_root.mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    archives = sorted(raw_dir.glob(args.pattern))
    if not archives:
        raise SystemExit(f"no archives found in {raw_dir} matching pattern '{args.pattern}'")

    ratio_tables: List[pd.DataFrame] = []
    for archive in archives:
        ratios = process_archive(archive, extract_root, force_extract=args.force_extract)
        if not ratios.empty:
            ratio_tables.append(ratios)

    if not ratio_tables:
        raise SystemExit("no redox ratios were generated; verify archive contents")

    combined = pd.concat(ratio_tables)
    ratio_columns = [col for col in combined.columns if col in RATIOS.keys()]
    ratio_columns = sorted(set(ratio_columns))
    combined = combined.loc[:, ratio_columns + ["dataset"]]
    combined.sort_index(inplace=True)
    combined.to_csv(output_path, sep="\t")

    logging.info("wrote redox ratios for %d samples to %s", combined.shape[0], output_path)


if __name__ == "__main__":
    main()
