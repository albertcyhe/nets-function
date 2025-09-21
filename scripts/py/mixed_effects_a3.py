#!/usr/bin/env python3
"""Mixed-effects models with organ interaction (Serpin -> endpoints).

Model per endpoint:
  value ~ Serpin_score_core * organ + log_total_PSMs + (1 | dataset)

Endpoints: footprint_index, THBS1_log_count, THBS1_cleave_idx, Proteo_DeltaFM.

Outputs results/tables/mixed_effects_a3.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--datasets", nargs="+", default=["PXD046330", "PXD005719", "PXD051579"])
    ap.add_argument("--processed-root", type=Path, default=Path("data/processed/proteomics"))
    ap.add_argument("--serpin", type=Path, default=Path("results/tables/serpin_scores.tsv"))
    ap.add_argument("--footprints", type=Path, default=Path("results/tables/footprints.tsv"))
    ap.add_argument("--thbs1cleave", type=Path, default=Path("results/tables/thbs1_cleave_idx.tsv"))
    ap.add_argument("--out", type=Path, default=Path("results/tables/mixed_effects_a3.tsv"))
    return ap.parse_args()


def derive_organ(dataset: str, sample: str) -> str:
    s = sample.upper()
    if dataset == "PXD005719":
        return "Brain" if "BR" in s else "NonBrain"
    if dataset == "PXD051579":
        return "NonBrain"
    if dataset == "PXD046330":
        return "NonBrain"
    return "NonBrain"


def load_total_psms(processed_root: Path, datasets: List[str]) -> pd.DataFrame:
    rows = []
    for ds in datasets:
        path = processed_root / ds / "gene_psm_matrix.tsv"
        if not path.exists():
            continue
        mat = pd.read_csv(path, sep="\t", index_col=0)
        counts = np.power(2.0, mat) - 1.0
        if counts.shape[0] < counts.shape[1] and {"SERPINA1","SERPINA3","SERPINI1"}.intersection(counts.columns):
            counts = counts.T
        total = counts.sum(axis=0)
        for s, val in total.items():
            rows.append({"dataset": ds, "sample": s, "log_total_PSMs": float(np.log1p(val))})
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    serpin = pd.read_csv(args.serpin, sep="\t")
    footprints = pd.read_csv(args.footprints, sep="\t")
    cleave = pd.read_csv(args.thbs1cleave, sep="\t") if args.thbs1cleave.exists() else pd.DataFrame()

    delta_rows = []
    for ds in args.datasets:
        dpath = args.processed_root / ds / "proteo_deltafm.tsv"
        if dpath.exists():
            d = pd.read_csv(dpath, sep="\t")[['Sample','Proteo_DeltaFM','THBS1_log_count']]
            d["dataset"] = ds
            d.rename(columns={"Sample": "sample"}, inplace=True)
            delta_rows.append(d)
    delta = pd.concat(delta_rows, ignore_index=True) if delta_rows else pd.DataFrame()

    covar = load_total_psms(args.processed_root, args.datasets)

    df = serpin.merge(footprints[['dataset','sample','footprint_index','control_footprint_index']], on=['dataset','sample'], how='left')
    df = df.merge(delta, on=['dataset','sample'], how='left')
    if not cleave.empty:
        if 'THBS1_cleave_idx' not in cleave.columns:
            cand = [c for c in cleave.columns if c.startswith('THBS1_cleave_idx')]
            if cand:
                cleave = cleave.rename(columns={cand[0]:'THBS1_cleave_idx'})
        df = df.merge(cleave[['dataset','sample','THBS1_cleave_idx']], on=['dataset','sample'], how='left')
    df = df.merge(covar, on=['dataset','sample'], how='left')
    df['organ'] = [derive_organ(d, s) for d, s in zip(df['dataset'], df['sample'])]

    endpoints = {
        'footprint_index': 'neg',
        'THBS1_log_count': 'pos',
        'THBS1_cleave_idx': 'neg',
        'Proteo_DeltaFM': 'neg'
    }

    rows = []
    for ycol, expect in endpoints.items():
        if ycol not in df.columns:
            continue
        data = df[['dataset','sample','Serpin_score_core','log_total_PSMs','organ',ycol]].replace([np.inf, -np.inf], np.nan).dropna()
        if data['organ'].nunique() < 2 or len(data) < 6:
            continue
        # Mixed-effects with dataset random intercept
        try:
            model = smf.mixedlm(f"{ycol} ~ Serpin_score_core * organ + log_total_PSMs", data, groups=data['dataset'])
            fit = model.fit(reml=False, method='lbfgs', maxiter=200)
            for term in ['Serpin_score_core','Serpin_score_core:organ[T.Brain]','log_total_PSMs','Intercept']:
                if term in fit.params.index:
                    rows.append({
                        'endpoint': ycol,
                        'expect': expect,
                        'term': term,
                        'estimate': float(fit.params[term]),
                        'se': float(fit.bse[term]),
                        'z': float(fit.tvalues[term]),
                        'p_value': float(fit.pvalues[term]),
                        'n': len(data)
                    })
        except Exception as exc:
            rows.append({'endpoint': ycol, 'term': 'model_error', 'estimate': np.nan, 'se': np.nan, 'z': np.nan, 'p_value': np.nan, 'n': len(data), 'error': str(exc)})

    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()
