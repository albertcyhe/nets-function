#!/usr/bin/env python3
"""Compute Proteo-ΔFM scores from Philosopher PSM reports.

The workflow expects a Philosopher workspace per dataset containing a
`psm.tsv` file. PSM-level evidence is collapsed to Gene × Sample spectral
counts (counts of high-confidence PSMs). Module scores are derived using the
proteomic NET-F and NET-M modules defined in `resources/modules/`.

Outputs (written to `data/processed/proteomics/`):

* `<dataset>/proteo_deltafm.tsv` – per-sample module scores and ΔFM values
* `<dataset>/gene_psm_matrix.tsv` – Gene × Sample matrix of log2(count+1)
* `proteo_deltafm_summary.tsv` – combined results across datasets
* `proteo_deltafm_vs_thbs1.tsv` – per-dataset correlation summary (Spearman)

Example
-------
```bash
micromamba run -n proteomics python scripts/py/proteomics_deltafm_score.py \
  --dataset PXD011796 --dataset PXD046330 --dataset PXD051579 \
  --workspace-root data/interim/proteomics \
  --output-root data/processed/proteomics \
  --module-f resources/modules/net_f_proteome.tsv \
  --module-m resources/modules/net_m_proteome.tsv
```
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats


PSM_FILE = "psm.tsv"
PROBABILITY_THRESHOLD = 0.9


@dataclass
class ModuleScores:
    dataset: str
    sample: str
    score_f: float
    score_m: float
    score_f_z: float
    score_m_z: float
    deltafm: float
    n_f_genes: int
    n_m_genes: int
    thbs1_count: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        action="append",
        dest="datasets",
        required=True,
        help="Dataset accession (repeat for multiple datasets)",
    )
    parser.add_argument(
        "--workspace-root",
        type=Path,
        default=Path("data/interim/proteomics"),
        help="Base directory containing <dataset>/fragger_closed/",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("data/processed/proteomics"),
        help="Base directory for Module outputs",
    )
    parser.add_argument(
        "--module-f",
        type=Path,
        default=Path("resources/modules/net_f_proteome.tsv"),
        help="TSV listing NET-F proteomic module genes (column 'Gene')",
    )
    parser.add_argument(
        "--module-m",
        type=Path,
        default=Path("resources/modules/net_m_proteome.tsv"),
        help="TSV listing NET-M proteomic module genes (column 'Gene')",
    )
    parser.add_argument(
        "--prob-threshold",
        type=float,
        default=PROBABILITY_THRESHOLD,
        help="Minimum PeptideProphet probability for PSM inclusion",
    )
    return parser.parse_args()


def load_module_genes(path: Path) -> List[str]:
    df = pd.read_csv(path, sep="\t")
    if "Gene" not in df.columns:
        raise ValueError(f"Module file {path} missing 'Gene' column")
    genes = df["Gene"].dropna().astype(str).str.strip().tolist()
    return sorted(set(genes))


def clean_sample_name(raw: str) -> str:
    name = Path(str(raw)).name
    if name.startswith("interact-"):
        name = name[len("interact-") :]
    if name.endswith(".pep.xml"):
        name = name[: -len(".pep.xml")]
    if name.endswith(".mzML"):
        name = name[: -len(".mzML")]
    return name


def derive_sample_column(df: pd.DataFrame) -> pd.Series:
    if "Spectrum File" in df.columns:
        unique_files = df["Spectrum File"].dropna().unique()
        if len(unique_files) > 1:
            return df["Spectrum File"].map(clean_sample_name)

    if "Spectrum" in df.columns:
        return df["Spectrum"].astype(str).str.split(".").str[0]

    raise ValueError("PSM table lacks both 'Spectrum File' diversity and 'Spectrum' column for sample inference")


def load_psm_counts(psm_path: Path, prob_threshold: float) -> pd.DataFrame:
    if not psm_path.exists():
        raise FileNotFoundError(f"PSM file not found: {psm_path}")

    df = pd.read_csv(psm_path, sep="\t")
    required_cols = {"Spectrum File", "Gene", "PeptideProphet Probability"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"PSM file {psm_path} missing columns: {sorted(missing)}")

    df = df[df["PeptideProphet Probability"] >= prob_threshold]
    df = df[df["Gene"].notna()]
    df["Sample"] = derive_sample_column(df)
    df["Gene"] = df["Gene"].astype(str).str.strip()
    df = df[df["Gene"] != ""]

    counts = df.groupby(["Sample", "Gene"]).size().unstack(fill_value=0)
    counts = counts.sort_index(axis=0)
    counts = counts.sort_index(axis=1)
    return counts


def zscore(series: pd.Series) -> pd.Series:
    if series.empty:
        return series
    if math.isclose(series.std(ddof=0), 0.0):
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (series - series.mean()) / series.std(ddof=0)


def compute_module_scores(
    dataset: str,
    counts: pd.DataFrame,
    genes_f: Sequence[str],
    genes_m: Sequence[str],
) -> Tuple[pd.DataFrame, List[ModuleScores]]:
    log_counts = np.log2(counts + 1)

    present_f = [g for g in genes_f if g in log_counts.columns]
    present_m = [g for g in genes_m if g in log_counts.columns]

    module_df = pd.DataFrame(index=log_counts.index)

    if present_f:
        module_df["score_F_raw"] = log_counts[present_f].mean(axis=1)
        module_df["n_F_genes"] = (log_counts[present_f] > 0).sum(axis=1)
    else:
        module_df["score_F_raw"] = np.nan
        module_df["n_F_genes"] = 0

    if present_m:
        module_df["score_M_raw"] = log_counts[present_m].mean(axis=1)
        module_df["n_M_genes"] = (log_counts[present_m] > 0).sum(axis=1)
    else:
        module_df["score_M_raw"] = np.nan
        module_df["n_M_genes"] = 0

    module_df["score_F_z"] = zscore(module_df["score_F_raw"].dropna()).reindex(
        module_df.index
    )
    module_df["score_M_z"] = zscore(module_df["score_M_raw"].dropna()).reindex(
        module_df.index
    )
    module_df["Proteo_DeltaFM"] = module_df["score_F_z"] - module_df["score_M_z"]

    thbs1 = log_counts.get("THBS1", pd.Series(0, index=log_counts.index))
    module_df["THBS1_log_count"] = thbs1

    records: List[ModuleScores] = []
    for sample, row in module_df.iterrows():
        records.append(
            ModuleScores(
                dataset=dataset,
                sample=sample,
                score_f=row["score_F_raw"],
                score_m=row["score_M_raw"],
                score_f_z=row["score_F_z"],
                score_m_z=row["score_M_z"],
                deltafm=row["Proteo_DeltaFM"],
                n_f_genes=int(row["n_F_genes"]),
                n_m_genes=int(row["n_M_genes"]),
                thbs1_count=row["THBS1_log_count"],
            )
        )

    return module_df, records


def spearman_summary(df: pd.DataFrame, dataset: str) -> Dict[str, float]:
    valid = df.dropna(subset=["Proteo_DeltaFM", "THBS1_log_count"])
    if len(valid) < 3:
        return {
            "dataset": dataset,
            "n": len(valid),
            "rho": np.nan,
            "p_value": np.nan,
        }
    rho, p_val = stats.spearmanr(
        valid["Proteo_DeltaFM"], valid["THBS1_log_count"]
    )
    return {
        "dataset": dataset,
        "n": int(len(valid)),
        "rho": float(rho),
        "p_value": float(p_val),
    }


def main() -> None:
    args = parse_args()
    genes_f = load_module_genes(args.module_f)
    genes_m = load_module_genes(args.module_m)

    summary_records: List[Dict[str, float]] = []
    combined_records: List[ModuleScores] = []

    for dataset in args.datasets:
        workspace = args.workspace_root / dataset / "fragger_closed"
        psm_path = workspace / PSM_FILE
        counts = load_psm_counts(psm_path, args.prob_threshold)

        module_df, records = compute_module_scores(dataset, counts, genes_f, genes_m)
        combined_records.extend(records)

        outdir = args.output_root / dataset
        outdir.mkdir(parents=True, exist_ok=True)
        counts.to_csv(outdir / "gene_psm_matrix.tsv", sep="\t")
        module_df.to_csv(outdir / "proteo_deltafm.tsv", sep="\t")

        summary_records.append(spearman_summary(module_df, dataset))

    if combined_records:
        combined_df = pd.DataFrame([r.__dict__ for r in combined_records])
        combined_df.to_csv(
            args.output_root / "proteo_deltafm_summary.tsv", sep="\t", index=False
        )

    if summary_records:
        summary_df = pd.DataFrame(summary_records)
        summary_df.to_csv(
            args.output_root / "proteo_deltafm_vs_thbs1.tsv", sep="\t", index=False
        )


if __name__ == "__main__":
    main()
