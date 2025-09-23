#!/usr/bin/env python3
"""Mixed-effects models for growth/bypass proxies.

Models (non-brain, errors handled):
  - Footprint_index ~ Serpin_score_core + log_total_PSMs + (1 | dataset)
  - THBS1_cleave_idx ~ Serpin_score_core + log_total_PSMs + (1 | dataset)
  - THBS1_log_count ~ Serpin_score_core + log_total_PSMs + (1 | dataset)
  - E2F_G2M_score ~ footprint_index + THBS1_cleave_idx + log_total_PSMs + (1 | dataset)

Brain: (n=1) -> store summary note.

Outputs results/tables/mixed_effects_growth.tsv and results/tables/mixed_effects_messages.txt
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
    ap.add_argument("--datasets", nargs="+", default=["PXD046330","PXD005719","PXD051579"])
    ap.add_argument("--processed-root", type=Path, default=Path("data/processed/proteomics"))
    ap.add_argument("--serpin", type=Path, default=Path("results/tables/serpin_scores.tsv"))
    ap.add_argument("--footprints", type=Path, default=Path("results/tables/footprints.tsv"))
    ap.add_argument("--thbs1cleave", type=Path, default=Path("results/tables/thbs1_cleave_idx.tsv"))
    ap.add_argument("--modules", type=Path, default=Path("results/tables/module_scores.tsv"))
    ap.add_argument("--out", type=Path, default=Path("results/tables/mixed_effects_growth.tsv"))
    ap.add_argument("--log", type=Path, default=Path("results/tables/mixed_effects_messages.txt"))
    return ap.parse_args()


def derive_organ(dataset: str, sample: str) -> str:
    s = sample.upper()
    if dataset == "PXD005719" and "BR" in s:
        return "Brain"
    if dataset in ("PXD046330","PXD051579"):
        return "NonBrain"
    return "NonBrain"


def load_total_psms(processed_root: Path, datasets: List[str]) -> pd.DataFrame:
    records = []
    for ds in datasets:
        path = processed_root / ds / 'gene_psm_matrix.tsv'
        if not path.exists():
            continue
        mat = pd.read_csv(path, sep='\t', index_col=0)
        counts = np.power(2.0, mat) - 1.0
        if counts.shape[0] < counts.shape[1] and {"SERPINA1","SERPINA3","SERPINI1"}.intersection(counts.columns):
            counts = counts.T
        total = counts.sum(axis=0)
        for sample, value in total.items():
            records.append({'dataset': ds, 'sample': sample, 'log_total_PSMs': float(np.log1p(value))})
    return pd.DataFrame(records)


def main() -> None:
    args = parse_args()
    serpin = pd.read_csv(args.serpin, sep='\t')
    footprints = pd.read_csv(args.footprints, sep='\t')
    modules = pd.read_csv(args.modules, sep='\t')
    cleave = pd.read_csv(args.thbs1cleave, sep='\t') if args.thbs1cleave.exists() else pd.DataFrame()

    delta_rows = []
    for ds in args.datasets:
        dpath = args.processed_root / ds / 'proteo_deltafm.tsv'
        if dpath.exists():
            d = pd.read_csv(dpath, sep='\t')[['Sample','Proteo_DeltaFM','THBS1_log_count']]
            d['dataset'] = ds
            d.rename(columns={'Sample':'sample'}, inplace=True)
            delta_rows.append(d)
    delta = pd.concat(delta_rows, ignore_index=True) if delta_rows else pd.DataFrame()

    covar = load_total_psms(args.processed_root, args.datasets)

    df = serpin.merge(footprints[['dataset','sample','footprint_index']], on=['dataset','sample'], how='left')
    if not cleave.empty:
        if 'THBS1_cleave_idx' not in cleave.columns:
            cand = [c for c in cleave.columns if c.startswith('THBS1_cleave_idx')]
            if cand:
                cleave = cleave.rename(columns={cand[0]:'THBS1_cleave_idx'})
        df = df.merge(cleave[['dataset','sample','THBS1_cleave_idx']], on=['dataset','sample'], how='left')
    df = df.merge(delta, on=['dataset','sample'], how='left')
    df = df.merge(modules[['dataset','sample','E2F_G2M_score','IL6_STAT3_score','CoOption_score']], on=['dataset','sample'], how='left')
    df = df.merge(covar, on=['dataset','sample'], how='left')
    df['organ'] = [derive_organ(d, s) for d, s in zip(df['dataset'], df['sample'])]

    models = [
        ('footprint_index', 'footprint_index ~ Serpin_score_core + log_total_PSMs + (1|dataset)'),
        ('THBS1_cleave_idx', 'THBS1_cleave_idx ~ Serpin_score_core + log_total_PSMs + (1|dataset)'),
        ('THBS1_log_count', 'THBS1_log_count ~ Serpin_score_core + log_total_PSMs + (1|dataset)'),
        ('E2F_G2M_score', 'E2F_G2M_score ~ footprint_index + THBS1_cleave_idx + log_total_PSMs + (1|dataset)')
    ]

    rows = []
    messages: List[str] = []

    for endpoint, formula in models:
        nonbrain = df[df['organ']=='NonBrain'][['dataset','sample','Serpin_score_core','log_total_PSMs','footprint_index','THBS1_cleave_idx','THBS1_log_count','Proteo_DeltaFM','E2F_G2M_score']]
        nonbrain = nonbrain.replace([np.inf,-np.inf], np.nan).dropna()
        if endpoint not in nonbrain.columns or nonbrain.empty:
            messages.append(f"{endpoint}: insufficient nonbrain data")
            continue
        target = endpoint
        exog_cols = ['Serpin_score_core','log_total_PSMs']
        if endpoint == 'E2F_G2M_score':
            exog_cols = ['footprint_index','THBS1_cleave_idx','log_total_PSMs']
        data = nonbrain[['dataset',target] + exog_cols].dropna().copy()
        data = data.rename(columns={target:'y'})
        for col in exog_cols:
            data.rename(columns={col:f"x_{col}"}, inplace=True)
        # ensure numeric
        for col in [c for c in data.columns if c.startswith('x_') or c=='y']:
            data[col] = pd.to_numeric(data[col], errors='coerce')
        data = data.dropna()
        if data.empty:
            messages.append(f"{endpoint}: no numeric data after coercion")
            continue
        formula_mixed = "y ~ " + " + ".join([f"x_{col}" for col in exog_cols])
        try:
            model = smf.mixedlm(formula_mixed, data, groups=data['dataset'])
            fit = model.fit(reml=False, method='lbfgs', maxiter=200)
            for col in exog_cols:
                name = f"x_{col}"
                if name in fit.params.index:
                    rows.append({
                        'endpoint': endpoint,
                        'term': col,
                        'estimate': float(fit.params[name]),
                        'se': float(fit.bse[name]),
                        'z': float(fit.tvalues[name]),
                        'p_value': float(fit.pvalues[name]),
                        'n': len(data)
                    })
        except Exception as exc:
            messages.append(f"{endpoint}: {exc}")

    if df[df['organ']=='Brain'].shape[0] < 3:
        messages.append("Brain stratum underpowered (n<3); only descriptive plots available.")

    outdf = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    outdf.to_csv(args.out, sep='\t', index=False)
    if messages:
        with open(args.log, 'w') as fh:
            for m in messages:
                fh.write(m + '\n')


if __name__ == '__main__':
    main()
