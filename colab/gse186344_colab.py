#!/usr/bin/env python3
"""
Colab runner for GSE186344 scRNA-seq pipeline (download via aria2c, preprocess, neutrophil low-ΔFM, THBS1 association).

Usage (in Colab notebook cell):

  !python colab/gse186344_colab.py \
      --gsm GSM5645895 GSM5645896 GSM5645897 GSM5645898 GSM5645899 GSM5645900 \
      GSM5645901 GSM5645902 GSM5645903 GSM5645904 GSM5645905 GSM5645906 \
      --project_dir "/content/drive/MyDrive/group B"

This script will:
  - apt-get install aria2 (Colab)
  - pip install a compatible Scanpy stack (numpy<2, scanpy 1.9.8, anndata 0.9.2, scikit-misc)
  - download each GSM’s 10x files (prefer h5; fallback to mtx) via aria2c in parallel
  - preprocess and save h5ad under <project_dir>/data/processed/scRNA/gse186344_colab.h5ad
  - compute neutrophil low-ΔFM subpopulation (NET-F/NET-M modules from project repo)
  - per-sample association: lowΔFM fraction vs tumor THBS1
  - write results under <project_dir>/results/tables
"""

import argparse
import os
import re
import sys
import subprocess
from pathlib import Path
from typing import List, Tuple

def run(cmd: List[str], check: bool=True, **kwargs):
    print("$", " ".join(cmd))
    r = subprocess.run(cmd, check=check, **kwargs)
    return r

def ensure_colab_env():
    # apt
    try:
        run(["apt-get","update"]) ; run(["apt-get","install","-y","aria2"])
    except Exception as e:
        print("[warn] apt-get failed or not available:", e)
    # pip stack (compatible)
    py = sys.executable
    pkgs = [
        "numpy<2", "pandas<2.2", "scanpy==1.9.8", "anndata==0.9.2", "scikit-misc",
        "matplotlib", "seaborn", "scikit-learn", "harmonypy"
    ]
    run([py,"-m","pip","install","-q","--upgrade"]+pkgs)

def fetch_suppl_for_gsm(gsm: str) -> List[str]:
    # GSM group prefix rule: GSM5645895 -> GSM5645nnn
    m = re.match(r"GSM(\d+)", gsm)
    assert m, f"Invalid GSM: {gsm}"
    prefix = m.group(1)[:4]
    grp = f"GSM{prefix}nnn"
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{grp}/{gsm}/suppl/"
    # Fetch listing
    import requests
    html = requests.get(base, timeout=60).text
    hrefs = re.findall(r'href="([^"]+)"', html)
    # Prefer single h5 (.h5 or .h5.gz)
    h5s = [h for h in hrefs if h.endswith('.h5') or h.endswith('.h5.gz')]
    urls = []
    if h5s:
        urls.append(base + h5s[0])
        return urls
    # Else search mtx triple
    needed = ['barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz']
    for f in needed:
        cand = [h for h in hrefs if h.endswith(f)]
        if cand:
            urls.append(base + cand[0])
    return urls

def download_gsms(gsms: List[str], outdir: Path) -> List[Tuple[str,Path]]:
    outdir.mkdir(parents=True, exist_ok=True)
    url_list = []
    for gsm in gsms:
        urls = fetch_suppl_for_gsm(gsm)
        if not urls:
            print(f"[warn] No files found for {gsm}")
            continue
        # put each GSM in its folder
        gdir = outdir / gsm
        gdir.mkdir(exist_ok=True)
        for u in urls:
            url_list.append((u,gdir))
    # Write a temp url list per target dir and fetch
    for gsm in gsms:
        pass
    # aria2c does not support per-line output dir; run per-GSM to keep simple
    for gsm in gsms:
        gdir = outdir / gsm
        # Collect urls for this gsm
        us = [u for u,d in url_list if d==gdir]
        if not us:
            continue
        tmp = gdir / "urls.txt"
        tmp.write_text("\n".join(us)+"\n")
        run(["aria2c","-c","-x","16","-s","16","--max-tries=0","--retry-wait=30","-d",str(gdir),"-i",str(tmp)])
    return [(gsm, outdir/gsm) for gsm in gsms]

def preprocess(samples: List[Tuple[str,Path]], project_dir: Path):
    import scanpy as sc, anndata as ad
    import numpy as np, pandas as pd
    from scipy import sparse as sp
    def read_sample(gsm: str, gdir: Path):
        # prefer h5 (.h5 or .h5.gz)
        h5s = list(gdir.glob('*.h5'))
        if not h5s:
            h5gz = list(gdir.glob('*.h5.gz'))
            if h5gz:
                import gzip, shutil
                src = h5gz[0]
                dst = gdir / src.name.replace('.h5.gz','.h5')
                with gzip.open(src, 'rb') as f_in, open(dst, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                h5s = [dst]
        if h5s:
            a = sc.read_10x_h5(str(h5s[0]))
            a.var_names = a.var_names.str.upper(); a.var_names_make_unique(); a.obs_names_make_unique()
            a.obs['sample']=gsm
            return a
        # else mtx (allow slight name variations)
        mtx_dir = gdir
        a = sc.read_10x_mtx(str(mtx_dir), var_names='gene_symbols')
        a.var_names = a.var_names.str.upper(); a.var_names_make_unique(); a.obs_names_make_unique()
        a.obs['sample']=gsm
        return a
    adatas=[]
    for gsm,gdir in samples:
        try:
            a = read_sample(gsm,gdir)
            adatas.append(a)
        except Exception as e:
            print(f"[warn] skip {gsm}: {e}")
    # Concatenate
    keys=[x.obs['sample'].iloc[0] for x in adatas]
    adata = ad.concat(adatas, join='outer', label='batch', keys=keys, fill_value=0)
    # QC
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs['pct_counts_mt'] < 20].copy()
    # Norm/log
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # HVG (fallbacks)
    n_top = 3000 if adata.n_vars>3500 else max(1000, int(adata.n_vars*0.8))
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor='seurat_v3', subset=True)
    except Exception:
        try:
            sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor='seurat', subset=True)
        except Exception:
            sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor='cell_ranger', subset=True)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)
    try:
        import harmonypy
        sc.external.pp.harmony_integrate(adata, key='sample')
        sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=20)
    except Exception:
        pass
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)
    out_h5ad = project_dir/"data/processed/scRNA/gse186344_colab.h5ad"
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(out_h5ad), compression='gzip')
    print("Wrote:", out_h5ad)
    return adata

def neutrophil_low_dfm_and_assoc(adata, project_dir: Path, neut_q: float = 0.70, epi_q: float = 0.75):
    import numpy as np, pandas as pd, scanpy as sc
    NEUT = ['S100A8','S100A9','CSF3R','FCGR3B','CXCR2','MPO','ELANE','LCN2','CEACAM8','MMP8','MMP9']
    T_NK = ['CD3D','CD3E','TRAC','NKG7','KLRD1']
    MONO = ['LYZ','S100A8','S100A9','CTSS','FCGR3A']
    TUM_EPI = ['EPCAM','KRT8','KRT18','KRT19','KRT7']
    # Score
    def score(adata, genes, key):
        genes=[g for g in genes if g in adata.var_names]
        sc.tl.score_genes(adata, gene_list=genes, score_name=key, use_raw=False)
    score(adata, NEUT, 'score_neut')
    score(adata, T_NK, 'score_tnk')
    score(adata, MONO, 'score_mono')
    # Gate neut: high neut, low T/NK, low Mono
    thr_neut = np.nanpercentile(adata.obs['score_neut'], neut_q*100)
    thr_tnk  = np.nanpercentile(adata.obs['score_tnk'].fillna(0), 60)
    thr_mono = np.nanpercentile(adata.obs['score_mono'].fillna(0), 60)
    mask = (adata.obs['score_neut'] >= thr_neut) & \
           (adata.obs['score_tnk'].fillna(0) <= thr_tnk) & \
           (adata.obs['score_mono'].fillna(0) <= thr_mono)
    neut = adata[mask].copy()
    # Modules
    net_f = [x.strip().upper() for x in (project_dir/"resources/modules/net_f_v1.tsv").read_text().splitlines() if x.strip()]
    net_m = [x.strip().upper() for x in (project_dir/"resources/modules/net_m_v2.tsv").read_text().splitlines() if x.strip()]
    score(neut, net_f, 'score_F') ; score(neut, net_m, 'score_M')
    z = lambda v: (v - np.nanmean(v))/(np.nanstd(v)+1e-8)
    neut.obs['zF']=z(neut.obs['score_F']); neut.obs['zM']=z(neut.obs['score_M'])
    neut.obs['deltaFM_cell']=neut.obs['zF']-neut.obs['zM']
    sc.pp.neighbors(neut, n_neighbors=15)
    sc.tl.leiden(neut, resolution=0.6, key_added='leiden_neut')
    cl_means = neut.obs.groupby('leiden_neut')['deltaFM_cell'].mean().sort_values()
    low_cl = cl_means.index[0]
    neut.obs['lowDFM_neut']=(neut.obs['leiden_neut']==low_cl).astype(int)
    frac = neut.obs.groupby('sample')['lowDFM_neut'].mean().rename('lowDFM_neut_fraction').reset_index()
    # Save fraction
    outdir = project_dir/"results/tables"; outdir.mkdir(parents=True, exist_ok=True)
    frac_path = outdir/"gse186344_neutrophil_lowDFM_fraction.tsv"
    frac.to_csv(frac_path, sep='\t', index=False); print("Wrote:", frac_path)
    # markers of low-ΔFM
    neut.obs['lowDFM_neut'] = neut.obs['lowDFM_neut'].astype('category')
    neut.obs['lowDFM_neut'] = neut.obs['lowDFM_neut'].cat.rename_categories({0:'neg',1:'pos'})
    sc.tl.rank_genes_groups(neut, 'lowDFM_neut', method='wilcoxon')
    markers = sc.get.rank_genes_groups_df(neut, group='pos')
    markers_path = outdir/"gse186344_lowDFM_markers.tsv"
    markers.to_csv(markers_path, sep='\t', index=False); print("Wrote:", markers_path)
    # Tumor epi and THBS1 association
    score(adata, TUM_EPI, 'score_epith')
    thr_e = np.nanpercentile(adata.obs['score_epith'], epi_q*100)
    tum = adata[adata.obs['score_epith']>=thr_e].copy()
    up = tum.var_names.str.upper(); thb_idx = np.where(up=='THBS1')[0]
    out_assoc = outdir/"gse186344_subpop_thbs1_assoc.tsv"
    if thb_idx.size==0:
        pd.DataFrame([{'rho':None,'p':None,'n':0,'join_key':'sample','note':'THBS1 not found'}]).to_csv(out_assoc, sep='\t', index=False)
        print("THBS1 not found; wrote NA association:", out_assoc)
    else:
        import numpy as np
        thb_name = tum.var_names[thb_idx[0]]
        X = tum[:,thb_name].X
        thb = X.toarray().flatten() if hasattr(X,'toarray') else (X.A1 if hasattr(X,'A1') else np.array(X).flatten())
        df = pd.DataFrame({'sample':tum.obs['sample'].values,'THBS1':thb})
        thbs1_per = df.groupby('sample')['THBS1'].median().rename('THBS1_median').reset_index()
        joined = frac.merge(thbs1_per, on='sample', how='inner')
        from scipy.stats import spearmanr
        rho,p = spearmanr(joined['lowDFM_neut_fraction'], joined['THBS1_median'], nan_policy='omit')
        pd.DataFrame([{'rho':rho,'p':p,'n':len(joined),'join_key':'sample','note':'Spearman lowΔFM vs THBS1'}]).to_csv(out_assoc, sep='\t', index=False)
        print("Wrote:", out_assoc)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gsm', nargs='+', required=True)
    ap.add_argument('--project_dir', default='/content/drive/MyDrive/group B')
    ap.add_argument('--neut_quantile', type=float, default=0.70, help='quantile for neutrophil gating (default 0.70)')
    ap.add_argument('--epi_quantile', type=float, default=0.75, help='quantile for epithelial (tumor-like) gating (default 0.75)')
    args = ap.parse_args()

    # Prepare env
    try:
        import google.colab  # type: ignore
        from google.colab import drive
        drive.mount('/content/drive')
    except Exception:
        print('[info] Not running in Colab or Drive already mounted.')
    # ensure_colab_env()  # disabled: you are using conda in Colab

    project_dir = Path(args.project_dir)
    raw_dir = Path('/content/gse186344_suppl')
    samples = download_gsms(args.gsm, raw_dir)
    adata = preprocess(samples, project_dir)
    neutrophil_low_dfm_and_assoc(adata, project_dir, neut_q=args.neut_quantile, epi_q=args.epi_quantile)

if __name__ == '__main__':
    main()
