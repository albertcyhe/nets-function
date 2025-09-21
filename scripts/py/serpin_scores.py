#!/usr/bin/env python3
"""Compute interpretable Serpin scores (core and extended) with CLR.

Core basket (secreted/extracellular): SERPINA1, SERPINA3, SERPINI1
Extended basket (sensitivity): + SLPI, ELAFIN (PI3), A2M

Normalization (label-free):
  - Start from gene_psm_matrix.tsv (log2(count+1)) as proxy if intensities absent
  - Convert to counts, add pseudocount 1, per-sample CLR: log(x) - mean(log(x))
  - Output two columns per dataset: Serpin_score_core, Serpin_score_ext

Outputs: results/tables/serpin_scores.tsv (dataset, sample, core, ext)
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


CORE = ["SERPINA1", "SERPINA3", "SERPINI1"]
EXT_ONLY = ["SLPI", "PI3", "ELAFIN", "A2M"]  # PI3 alias considered


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", action="append", dest="datasets", required=True)
    ap.add_argument("--processed-root", type=Path, default=Path("data/processed/proteomics"))
    ap.add_argument("--out", type=Path, default=Path("results/tables/serpin_scores.tsv"))
    return ap.parse_args()


def load_counts(ds_dir: Path) -> pd.DataFrame:
    mat = pd.read_csv(ds_dir / "gene_psm_matrix.tsv", sep="\t", index_col=0)
    counts = np.power(2.0, mat) - 1.0
    return counts


def clr(df: pd.DataFrame) -> pd.DataFrame:
    X = df.copy()
    X = X + 1.0
    logX = np.log(X)
    gm = logX.mean(axis=0)
    return logX.subtract(gm, axis=1)


def main() -> None:
    args = parse_args()
    rows = []
    for ds in args.datasets:
        ds_dir = args.processed_root / ds
        counts = load_counts(ds_dir)
        # Orient genes as rows
        if len(set(CORE + EXT_ONLY).intersection(counts.columns)) > len(set(CORE + EXT_ONLY).intersection(counts.index)):
            counts = counts.T
        z = clr(counts)
        core_genes = [g for g in CORE if g in z.index]
        ext_genes = core_genes + [g for g in EXT_ONLY if g in z.index]
        core = z.loc[core_genes].mean(axis=0) if core_genes else pd.Series(np.nan, index=z.columns)
        ext = z.loc[ext_genes].mean(axis=0) if ext_genes else pd.Series(np.nan, index=z.columns)
        for sample in z.columns:
            rows.append({
                "dataset": ds,
                "sample": sample,
                "Serpin_score_core": float(core.get(sample, np.nan)),
                "Serpin_score_ext": float(ext.get(sample, np.nan)),
            })
    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()

