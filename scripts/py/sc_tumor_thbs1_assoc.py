#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

EPITH_MARKERS = ['EPCAM','KRT8','KRT18','KRT19','KRT7']


def score(adata, genes, key):
    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        adata.obs[key] = np.nan
        return
    sc.tl.score_genes(adata, gene_list=genes, score_name=key, use_raw=False)


def main():
    ap = argparse.ArgumentParser(description='Associate neutrophil low-ΔFM fraction with tumor-cell THBS1 per sample.')
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--frac', required=True, help='Output from sc_neutrophil_subsets.py (per-sample fraction)')
    ap.add_argument('--out', default='results/tables/gse186344_subpop_thbs1_assoc.tsv')
    ap.add_argument('--epi_quantile', type=float, default=0.75, help='Quantile for epithelial (tumor-like) gating (default 0.75)')
    args = ap.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    adata.var_names = adata.var_names.str.upper()

    # identify tumor-like epithelial cells
    score(adata, EPITH_MARKERS, 'score_epith')
    thr = np.nanpercentile(adata.obs['score_epith'], args.epi_quantile*100)
    tum = adata[adata.obs['score_epith'] >= thr].copy()
    if tum.n_obs == 0:
        raise SystemExit('No tumor-like cells detected by epithelial markers')

    if 'sample' in tum.obs.columns:
        skey = 'sample'
    elif 'batch' in tum.obs.columns:
        skey = 'batch'
    else:
        raise SystemExit('No sample or batch key in obs')

    # per-sample THBS1 expression (log-normalized)
    # robust THBS1 locator (case-insensitive)
    thb = None
    up = tum.var_names.str.upper()
    if 'THBS1' in set(up):
        thb = tum.var_names[up.get_loc('THBS1')]
    if thb is None:
        # write NA association and exit gracefully
        Path(Path(args.out).parent).mkdir(parents=True, exist_ok=True)
        pd.DataFrame([{'rho': np.nan, 'p': np.nan, 'n': 0, 'join_key': skey, 'note': 'THBS1 not found; association NA'}]).to_csv(args.out, sep='\t', index=False)
        print('THBS1 not present; wrote NA association to', args.out)
        return
    thbs = tum[:, thb].X
    thbs_arr = thbs.toarray().flatten() if hasattr(thbs, 'toarray') else (thbs.A1 if hasattr(thbs, 'A1') else np.array(thbs).flatten())
    thbs1_per_sample = pd.DataFrame({skey: tum.obs[skey].values, 'THBS1': thbs_arr}).groupby(skey)['THBS1'].median().rename('THBS1_median').reset_index()

    frac = pd.read_csv(args.frac, sep='\t')
    join_key = skey
    df = pd.merge(frac, thbs1_per_sample, left_on=join_key, right_on=join_key, how='inner')
    if df.empty:
        raise SystemExit('No overlapping samples between fraction table and tumor THBS1 table')

    rho, p = spearmanr(df['lowDFM_neut_fraction'], df['THBS1_median'], nan_policy='omit')
    out = pd.DataFrame([{
        'rho': rho, 'p': p, 'n': df.shape[0],
        'join_key': join_key,
        'note': 'Spearman: lowDFM_neut_fraction vs tumor THBS1 median'
    }])
    Path(Path(args.out).parent).mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep='\t', index=False)
    print('Wrote:', args.out)


if __name__ == '__main__':
    main()
