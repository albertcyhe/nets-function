#!/usr/bin/env python3
"""Compute ISI (Inhibitor–Protease Stoichiometry Index) and associate with Proteo-ΔFM.

Inputs (per dataset):
  - data/processed/proteomics/<dataset>/gene_psm_matrix.tsv (log2(count+1))
  - data/processed/proteomics/<dataset>/proteo_deltafm.tsv
  - Optionally data/processed/proteomics/<dataset>/protein_annot.tsv (linear intensities)

Outputs:
  - results/tables/isi_per_sample.tsv – per-sample inhibitor/target sums and ISI
  - results/tables/isi_models.tsv – per-dataset OLS associations (ΔFM ~ ISI)
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm


DEFAULT_INHIBITORS_CORE = [
    "SERPINA1", "SERPINA3", "SERPINB1", "SERPINB6", "SERPINB8",
]
DEFAULT_INHIBITORS_EXT = [
    "SLPI", "PI3", "PI3/Elafin", "ELAFIN", "WFDC14", "A2M",
]
DEFAULT_TARGETS = ["ELANE", "PRTN3", "CTSG"]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", action="append", dest="datasets", required=True,
                    help="Dataset(s): e.g., PXD005719 (repeatable)")
    ap.add_argument("--processed-root", type=Path, default=Path("data/processed/proteomics"))
    ap.add_argument("--outdir", type=Path, default=Path("results/tables"))
    ap.add_argument("--inhibitors", type=Path, default=None,
                    help="TSV with column 'Gene' listing inhibitors (defaults included)")
    ap.add_argument("--targets", type=Path, default=None,
                    help="TSV with column 'Gene' listing protease targets (defaults included)")
    ap.add_argument("--include-extended", action="store_true",
                    help="Include extended inhibitors (SLPI/ELAFIN/A2M) for sensitivity")
    ap.add_argument("--prefer-protein-intensity", action="store_true",
                    help="Use protein_annot.tsv intensities if available (else fallback to gene_psm counts)")
    return ap.parse_args()


def load_gene_psm_matrix(ds_dir: Path) -> pd.DataFrame:
    path = ds_dir / "gene_psm_matrix.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Missing gene_psm_matrix: {path}")
    mat = pd.read_csv(path, sep="\t", index_col=0)
    # convert log2(count+1) back to counts
    counts = np.power(2.0, mat) - 1.0
    return counts


def load_protein_annot(ds_dir: Path) -> pd.DataFrame | None:
    path = ds_dir / "protein_annot.tsv"
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    # intensity columns by Philosopher naming
    inten_cols = [c for c in df.columns if c.endswith("Intensity")]
    if not inten_cols:
        return None
    # collapse to gene by sum of intensities across peptides
    pivot = df.groupby("Gene")[inten_cols].sum()
    # rename sample columns as sample ids
    # Philosopher report uses generic columns; treat as three columns if present
    pivot.columns = [
        c.replace(" Total Intensity", "").replace(" Unique Intensity", "").replace(" Razor Intensity", "")
        if " Intensity" in c else c for c in pivot.columns
    ]
    return pivot


def load_deltafm(ds_dir: Path) -> pd.DataFrame:
    path = ds_dir / "proteo_deltafm.tsv"
    df = pd.read_csv(path, sep="\t")
    return df


def get_gene_lists(args: argparse.Namespace) -> Tuple[List[str], List[str]]:
    inhibitors = list(DEFAULT_INHIBITORS_CORE)
    if args.include_extended:
        inhibitors = inhibitors + DEFAULT_INHIBITORS_EXT
    if args.inhibitors and args.inhibitors.exists():
        custom = pd.read_csv(args.inhibitors, sep="\t")
        if "Gene" in custom.columns:
            inhibitors = [g for g in custom["Gene"].astype(str).str.strip().tolist() if g]
    targets = list(DEFAULT_TARGETS)
    if args.targets and args.targets.exists():
        custom_t = pd.read_csv(args.targets, sep="\t")
        if "Gene" in custom_t.columns:
            targets = [g for g in custom_t["Gene"].astype(str).str.strip().tolist() if g]
    return inhibitors, targets


def derive_organ(dataset: str, sample: str) -> str:
    s = sample.upper()
    if dataset == "PXD005719":
        return "Brain" if "BR" in s else "Other"
    if dataset == "PXD046330":
        return "InVitro"
    if dataset == "PXD051579":
        return "InVitro"
    return "Unknown"


def main() -> None:
    args = parse_args()
    inhibitors, targets = get_gene_lists(args)
    rows_per_sample: List[Dict[str, object]] = []
    model_rows: List[Dict[str, object]] = []

    for ds in args.datasets:
        ds_dir = args.processed_root / ds
        # choose source matrix
        expr = None
        if args.prefer_protein_intensity:
            expr = load_protein_annot(ds_dir)
        if expr is None:
            expr = load_gene_psm_matrix(ds_dir)

        # Orient to gene x sample if necessary
        probe_genes = set(inhibitors + targets)
        n_match_index = len(probe_genes.intersection(expr.index))
        n_match_cols = len(probe_genes.intersection(expr.columns))
        if n_match_cols > n_match_index:
            expr = expr.T

        # intersect genes
        inh = [g for g in inhibitors if g in expr.index]
        tgt = [g for g in targets if g in expr.index]
        if not inh:
            inh = []
        if not tgt:
            tgt = []

        # sums per sample
        inh_sum = expr.loc[inh].sum(axis=0) if inh else pd.Series(0.0, index=expr.columns)
        tgt_sum = expr.loc[tgt].sum(axis=0) if tgt else pd.Series(0.0, index=expr.columns)
        isi = np.log2((inh_sum.replace(0, np.nan) + 1e-6) / (tgt_sum.replace(0, np.nan) + 1e-6))

        # join with ΔFM table
        delta = load_deltafm(ds_dir)
        delta = delta.set_index("Sample")

        for sample in expr.columns:
            d = {
                "dataset": ds,
                "sample": sample,
                "inhibitor_sum": float(inh_sum.get(sample, 0.0)),
                "target_sum": float(tgt_sum.get(sample, 0.0)),
                "ISI": float(isi.get(sample, np.nan)),
                "Proteo_DeltaFM": float(delta.loc[sample, "Proteo_DeltaFM"]) if sample in delta.index else np.nan,
                "organ": derive_organ(ds, sample),
            }
            rows_per_sample.append(d)

        # simple OLS per dataset (ΔFM ~ ISI), robust se
        df = pd.DataFrame([r for r in rows_per_sample if r["dataset"] == ds]).dropna(subset=["ISI", "Proteo_DeltaFM"])
        if len(df) >= 3:
            X = sm.add_constant(df[["ISI"]])
            y = df["Proteo_DeltaFM"]
            ols = sm.OLS(y, X).fit(cov_type="HC3")
            model_rows.append({
                "dataset": ds,
                "n": int(len(df)),
                "beta_ISI": float(ols.params.get("ISI", np.nan)),
                "beta_SE": float(ols.bse.get("ISI", np.nan)),
                "p_value": float(ols.pvalues.get("ISI", np.nan)),
                "r2": float(ols.rsquared),
            })

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows_per_sample).to_csv(outdir / "isi_per_sample.tsv", sep="\t", index=False)
    pd.DataFrame(model_rows).to_csv(outdir / "isi_models.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
