#!/usr/bin/env python3
"""Correlate module scores with key endpoints (brain vs non-brain).

Outputs results/tables/module_assoc.tsv with Spearman rho and p.

Endpoints considered:
  - Footprint_index
  - THBS1_cleave_idx
  - THBS1_log_count
  - Proteo_DeltaFM

Modules: E2F_G2M_score, IL6_STAT3_score, CoOption_score

Strata: Brain, NonBrain, Combined.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy import stats


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--module", type=Path, default=Path("results/tables/module_scores.tsv"))
    ap.add_argument("--serpin", type=Path, default=Path("results/tables/serpin_scores.tsv"))
    ap.add_argument("--footprints", type=Path, default=Path("results/tables/footprints.tsv"))
    ap.add_argument("--thbs1cleave", type=Path, default=Path("results/tables/thbs1_cleave_idx.tsv"))
    ap.add_argument("--processed-root", type=Path, default=Path("data/processed/proteomics"))
    ap.add_argument("--datasets", nargs="+", default=["PXD046330","PXD005719","PXD051579"])
    ap.add_argument("--out", type=Path, default=Path("results/tables/module_assoc.tsv"))
    return ap.parse_args()


def derive_organ(dataset: str, sample: str) -> str:
    s = sample.upper()
    if dataset == "PXD005719":
        return "Brain" if "BR" in s else "NonBrain"
    if dataset in ("PXD046330","PXD051579"):
        return "NonBrain"
    return "NonBrain"


def load_delta(processed_root: Path, datasets: List[str]) -> pd.DataFrame:
    frames = []
    for ds in datasets:
        path = processed_root / ds / "proteo_deltafm.tsv"
        if not path.exists():
            continue
        df = pd.read_csv(path, sep='\t')[['Sample','Proteo_DeltaFM','THBS1_log_count']]
        df['dataset'] = ds
        df.rename(columns={'Sample':'sample'}, inplace=True)
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main() -> None:
    args = parse_args()
    module = pd.read_csv(args.module, sep='\t')
    serpin = pd.read_csv(args.serpin, sep='\t')
    footprints = pd.read_csv(args.footprints, sep='\t')
    cleave = pd.read_csv(args.thbs1cleave, sep='\t') if args.thbs1cleave.exists() else pd.DataFrame()
    delta = load_delta(args.processed_root, args.datasets)

    df = module.merge(serpin[['dataset','sample','Serpin_score_core']], on=['dataset','sample'], how='left')
    df = df.merge(footprints[['dataset','sample','footprint_index']], on=['dataset','sample'], how='left')
    if not cleave.empty:
        if 'THBS1_cleave_idx' not in cleave.columns:
            cand = [c for c in cleave.columns if c.startswith('THBS1_cleave_idx')]
            if cand:
                cleave = cleave.rename(columns={cand[0]:'THBS1_cleave_idx'})
        df = df.merge(cleave[['dataset','sample','THBS1_cleave_idx']], on=['dataset','sample'], how='left')
    df = df.merge(delta, on=['dataset','sample'], how='left')
    df['organ'] = [derive_organ(d, s) for d, s in zip(df['dataset'], df['sample'])]

    modules = [c for c in df.columns if c.endswith('_score')]
    endpoints = ['footprint_index','THBS1_cleave_idx','THBS1_log_count','Proteo_DeltaFM']
    strata = ['Brain','NonBrain','Combined']
    rows = []
    for mod in modules:
        for y in endpoints:
            if y not in df.columns:
                continue
            for strat in strata:
                if strat == 'Combined':
                    sub = df[[mod,y,'Serpin_score_core','organ']].copy()
                else:
                    sub = df[df['organ']==strat][[mod,y,'Serpin_score_core','organ']]
                sub = sub.dropna()
                if len(sub) < 3:
                    rows.append({'module':mod,'endpoint':y,'stratum':strat,'rho':np.nan,'p_value':np.nan,'n':len(sub)})
                    continue
                rho, p = stats.spearmanr(sub[mod], sub[y])
                rows.append({'module':mod,'endpoint':y,'stratum':strat,'rho':rho,'p_value':p,'n':len(sub)})

    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()

