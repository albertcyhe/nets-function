#!/usr/bin/env python3
"""Compute CLR-based module scores (proliferation, IL6/STAT3, co-option).

Modules (hand curated, amend as needed):
  - E2F_G2M (proliferation): MKI67, PCNA, CCNB1, CDK1, MCM2-7, TOP2A, BIRC5
  - IL6_STAT3: IL6, STAT3, SOCS3, JAK1, JAK2, IL6ST, CEBPB, FOS
  - CoOption (vascular/astrocytic): GFAP, ITGA6, ITGB1, LAMC1, COL4A1, VEGFA, ANGPT1, ANGPT2, S100B

Outputs results/tables/module_scores.tsv with CLR means per module.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


MODULES = {
    "E2F_G2M_score": [
        "MKI67","PCNA","CCNB1","CDK1",
        "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7",
        "TOP2A","BIRC5"
    ],
    "IL6_STAT3_score": [
        "IL6","STAT3","SOCS3","JAK1","JAK2","IL6ST","CEBPB","FOS"
    ],
    "CoOption_score": [
        "GFAP","ITGA6","ITGB1","LAMC1","COL4A1","VEGFA","ANGPT1","ANGPT2","S100B"
    ]
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--datasets", nargs="+", default=["PXD046330","PXD005719","PXD051579"])
    ap.add_argument("--processed-root", type=Path, default=Path("data/processed/proteomics"))
    ap.add_argument("--out", type=Path, default=Path("results/tables/module_scores.tsv"))
    return ap.parse_args()


def load_counts(path: Path) -> pd.DataFrame:
    mat = pd.read_csv(path, sep='\t', index_col=0)
    counts = np.power(2.0, mat) - 1.0
    return counts


def clr(df: pd.DataFrame) -> pd.DataFrame:
    logX = np.log(df + 1.0)
    gm = logX.mean(axis=0)
    return logX.subtract(gm, axis=1)


def main() -> None:
    args = parse_args()
    rows = []
    for ds in args.datasets:
        path = args.processed_root / ds / 'gene_psm_matrix.tsv'
        if not path.exists():
            continue
        counts = load_counts(path)
        all_genes = set().union(*MODULES.values())
        if len(all_genes.intersection(counts.columns)) > len(all_genes.intersection(counts.index)):
            counts = counts.T
        Z = clr(counts)
        for sample in Z.columns:
            row = {'dataset': ds, 'sample': sample}
            for mod_name, genes in MODULES.items():
                present = [g for g in genes if g in Z.index]
                row[mod_name] = float(Z.loc[present, sample].mean()) if present else np.nan
            rows.append(row)
    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()

