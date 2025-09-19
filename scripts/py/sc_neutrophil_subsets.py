#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc


NEUT_MARKERS = [
    'S100A8','S100A9','CSF3R','FCGR3B','CXCR2','MPO','ELANE','LCN2','CEACAM8','MMP8','MMP9'
]
T_NK_MARKERS = ['CD3D','CD3E','TRAC','NKG7','KLRD1']
MONO_MARKERS = ['LYZ','S100A8','S100A9','CTSS','FCGR3A']


def score(adata, genes, key):
    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        adata.obs[key] = np.nan
        return
    sc.tl.score_genes(adata, gene_list=genes, score_name=key, use_raw=False)


def read_modules(fpath):
    g = [x.strip().upper() for x in Path(fpath).read_text().splitlines() if x.strip()]
    return g


def main():
    ap = argparse.ArgumentParser(description='Identify neutrophil low-ΔFM subpopulation and compute per-sample fractions.')
    ap.add_argument('--h5ad', required=True, help='Preprocessed h5ad from sc_preprocess_gse186344.py')
    ap.add_argument('--net_f', default='resources/modules/net_f_v1.tsv')
    ap.add_argument('--net_m', default='resources/modules/net_m_v2.tsv')
    ap.add_argument('--outdir', default='results/tables')
    ap.add_argument('--neut_quantile', type=float, default=0.70, help='Quantile for neutrophil gating (default 0.70)')
    args = ap.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    adata.var_names = adata.var_names.str.upper()

    # score lineage markers
    score(adata, NEUT_MARKERS, 'score_neut')
    score(adata, T_NK_MARKERS, 'score_tnk')
    score(adata, MONO_MARKERS, 'score_mono')

    # neutrophil gate: high neut (>= neut_quantile), and low T/NK & Mono (<= 60th percentile)
    thr_neut = np.nanpercentile(adata.obs['score_neut'], args.neut_quantile*100)
    thr_tnk = np.nanpercentile(adata.obs['score_tnk'].fillna(0), 60)
    thr_mono = np.nanpercentile(adata.obs['score_mono'].fillna(0), 60)
    mask = (adata.obs['score_neut'] >= thr_neut) & (adata.obs['score_tnk'].fillna(0) <= thr_tnk) & (adata.obs['score_mono'].fillna(0) <= thr_mono)
    neut = adata[mask].copy()
    # compute module scores in neutrophils
    net_f = read_modules(args.net_f)
    net_m = read_modules(args.net_m)
    score(neut, net_f, 'score_F')
    score(neut, net_m, 'score_M')
    # z per cell-population
    def z(v):
        v = np.asarray(v)
        return (v - np.nanmean(v)) / (np.nanstd(v) + 1e-8)
    neut.obs['zF'] = z(neut.obs['score_F'])
    neut.obs['zM'] = z(neut.obs['score_M'])
    neut.obs['deltaFM_cell'] = neut.obs['zF'] - neut.obs['zM']

    # cluster within neutrophils and identify lowest ΔFM cluster
    sc.pp.neighbors(neut, n_neighbors=15, use_rep='X_pca' if 'X_pca' in neut.obsm else None)
    sc.tl.leiden(neut, resolution=0.6, key_added='leiden_neut')
    cl_means = neut.obs.groupby('leiden_neut')['deltaFM_cell'].mean().sort_values()
    low_cl = cl_means.index[0]
    neut.obs['lowDFM_neut'] = (neut.obs['leiden_neut'] == low_cl).astype(int)
    neut.obs['lowDFM_neut_int'] = neut.obs['lowDFM_neut'].astype(int)
    neut.obs['lowDFM_neut'] = neut.obs['lowDFM_neut'].astype('category')
    # rename categories to strings for stable downstream selection
    if hasattr(neut.obs['lowDFM_neut'].cat, 'rename_categories'):
        neut.obs['lowDFM_neut'] = neut.obs['lowDFM_neut'].cat.rename_categories(lambda x: 'pos' if int(x)==1 else 'neg')

    # per-sample fractions
    if 'sample' not in neut.obs.columns:
        # fall back to batch if present
        sample_key = 'batch' if 'batch' in neut.obs.columns else None
    else:
        sample_key = 'sample'
    if sample_key is None:
        raise SystemExit('No sample or batch key in obs')
    frac = neut.obs.groupby(sample_key)['lowDFM_neut_int'].mean().rename('lowDFM_neut_fraction').reset_index()

    # marker genes for the low-ΔFM cluster
    sc.tl.rank_genes_groups(neut, 'lowDFM_neut', method='wilcoxon')
    markers = sc.get.rank_genes_groups_df(neut, group='pos')

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    frac_path = outdir / 'gse186344_neutrophil_lowDFM_fraction.tsv'
    frac.to_csv(frac_path, sep='\t', index=False)
    markers_path = outdir / 'gse186344_lowDFM_markers.tsv'
    markers.to_csv(markers_path, sep='\t', index=False)
    print(f'Wrote: {frac_path}')
    print(f'Wrote: {markers_path}')


if __name__ == '__main__':
    main()
