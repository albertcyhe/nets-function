#!/usr/bin/env python3
"""Compute negative control markers for association sanity checks.

Outputs per sample:
  - ECM_nc_score (CLR mean of laminin subunits): LAMA1/2/3, LAMB1/2/3, LAMC1/2/3
  - control_footprint_index: from footprints.tsv

Writes results/tables/neg_controls.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


LAMININS = [
    "LAMA1","LAMA2","LAMA3","LAMB1","LAMB2","LAMB3","LAMC1","LAMC2","LAMC3"
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--datasets", nargs="+", default=["PXD046330","PXD005719","PXD051579"])
    ap.add_argument("--processed-root", type=Path, default=Path("data/processed/proteomics"))
    ap.add_argument("--footprints", type=Path, default=Path("results/tables/footprints.tsv"))
    ap.add_argument("--out", type=Path, default=Path("results/tables/neg_controls.tsv"))
    return ap.parse_args()


def load_counts(path: Path) -> pd.DataFrame:
    mat = pd.read_csv(path, sep='\t', index_col=0)
    counts = np.power(2.0, mat) - 1.0
    return counts


def clr(df: pd.DataFrame) -> pd.DataFrame:
    X = df + 1.0
    logX = np.log(X)
    gm = logX.mean(axis=0)
    return logX.subtract(gm, axis=1)


def main() -> None:
    args = parse_args()
    foot = pd.read_csv(args.footprints, sep='\t')
    rows = []
    for ds in args.datasets:
        path = args.processed_root / ds / 'gene_psm_matrix.tsv'
        if not path.exists():
            continue
        counts = load_counts(path)
        # orient genes x samples
        if len(set(LAMININS).intersection(counts.columns)) > len(set(LAMININS).intersection(counts.index)):
            counts = counts.T
        Z = clr(counts)
        lam = [g for g in LAMININS if g in Z.index]
        ecm = Z.loc[lam].mean(axis=0) if lam else pd.Series(np.nan, index=Z.columns)
        for s in Z.columns:
            cfi = foot.loc[(foot['dataset']==ds) & (foot['sample']==s), 'control_footprint_index']
            rows.append({
                'dataset': ds,
                'sample': s,
                'ECM_nc_score': float(ecm.get(s, np.nan)),
                'control_footprint_index': float(cfi.iloc[0]) if len(cfi)>0 else np.nan,
            })
    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()

