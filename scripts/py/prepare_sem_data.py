#!/usr/bin/env python3
"""Assemble SEM input table from existing summaries."""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--serpin", type=Path, default=Path("results/tables/serpin_scores.tsv"))
    ap.add_argument("--footprints", type=Path, default=Path("results/tables/footprints.tsv"))
    ap.add_argument("--thbs1cleave", type=Path, default=Path("results/tables/thbs1_cleave_idx.tsv"))
    ap.add_argument("--processed-root", type=Path, default=Path("data/processed/proteomics"))
    ap.add_argument("--module", type=Path, default=Path("results/tables/module_scores.tsv"))
    ap.add_argument("--out", type=Path, default=Path("results/tables/sem_input.tsv"))
    return ap.parse_args()


def derive_organ(dataset: str, sample: str) -> str:
    s = sample.upper()
    if dataset == "PXD005719" and "BR" in s:
        return "Brain"
    if dataset in ("PXD046330", "PXD051579"):
        return "NonBrain"
    return "NonBrain"


def main() -> None:
    args = parse_args()
    serpin = pd.read_csv(args.serpin, sep='\t')
    footprints = pd.read_csv(args.footprints, sep='\t')
    cleave = pd.read_csv(args.thbs1cleave, sep='\t') if args.thbs1cleave.exists() else pd.DataFrame()
    modules = pd.read_csv(args.module, sep='\t') if args.module.exists() else pd.DataFrame()

    # proteo delta
    delta_frames = []
    for ds in serpin['dataset'].unique():
        path = args.processed_root / ds / 'proteo_deltafm.tsv'
        if path.exists():
            d = pd.read_csv(path, sep='\t')
            cols = ['Sample','Proteo_DeltaFM']
            if 'THBS1_log_count' in d.columns:
                cols.append('THBS1_log_count')
            d = d[cols].rename(columns={'Sample':'sample'})
            d['dataset'] = ds
            delta_frames.append(d)
    delta = pd.concat(delta_frames, ignore_index=True) if delta_frames else pd.DataFrame()

    df = serpin.merge(footprints[['dataset','sample','footprint_index']], on=['dataset','sample'], how='left')
    if not cleave.empty:
        if 'THBS1_cleave_idx' not in cleave.columns:
            cand = [c for c in cleave.columns if c.startswith('THBS1_cleave_idx')]
            if cand:
                cleave = cleave.rename(columns={cand[0]:'THBS1_cleave_idx'})
        df = df.merge(cleave[['dataset','sample','THBS1_cleave_idx']], on=['dataset','sample'], how='left')
    df = df.merge(delta, on=['dataset','sample'], how='left')
    if not modules.empty:
        df = df.merge(modules[['dataset','sample','E2F_G2M_score','IL6_STAT3_score','CoOption_score']], on=['dataset','sample'], how='left')

    df['organ'] = [derive_organ(d, s) for d, s in zip(df['dataset'], df['sample'])]
    df = df.dropna(subset=['Serpin_score_core','footprint_index','THBS1_cleave_idx','Proteo_DeltaFM'])
    df.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()

