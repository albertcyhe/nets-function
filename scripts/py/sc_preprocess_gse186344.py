#!/usr/bin/env python3
import argparse
import os
import tarfile
import gzip
import tempfile
from pathlib import Path
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp


def find_10x_dirs(root: Path):
    tenx = []
    for p in root.rglob('*'):
        if p.is_dir():
            mtx = list(p.glob('matrix.mtx*'))
            feats = list(p.glob('features.tsv*')) + list(p.glob('genes.tsv*'))
            bar = list(p.glob('barcodes.tsv*'))
            if mtx and feats and bar:
                tenx.append(p)
    return sorted(set(tenx))


def find_10x_h5(root: Path):
    h5s = []
    for p in root.rglob('*.h5'):
        if 'feature_bc_matrix' in p.name.lower() or 'filtered_feature_bc_matrix' in p.as_posix().lower():
            h5s.append(p)
    return sorted(set(h5s))


def read_10x_dir(p: Path):
    try:
        adata = sc.read_10x_mtx(p.as_posix(), var_names='gene_symbols', cache=True)
    except Exception:
        adata = sc.read_10x_mtx(p.as_posix(), var_names='gene_ids', cache=True)
    # ensure gene symbols uppercase where available
    if 'gene_symbols' in adata.var.columns:
        adata.var_names = adata.var['gene_symbols']
    adata.var_names = adata.var_names.str.upper()
    adata.var_names_make_unique()
    # ensure unique cell barcodes
    adata.obs_names_make_unique()
    # annotate sample from directory name
    adata.obs['sample'] = p.parts[-2] if p.name.lower() in {'filtered_feature_bc_matrix','outs'} else p.name
    return adata


def read_10x_h5(fp: Path):
    adata = sc.read_10x_h5(fp.as_posix())
    # 10x h5 reader usually sets symbols; enforce uppercase
    adata.var_names = adata.var_names.str.upper()
    adata.var_names_make_unique()
    # ensure unique cell barcodes
    adata.obs_names_make_unique()
    # sample inferred from parent directory or stem
    adata.obs['sample'] = fp.parent.name or fp.stem
    return adata


def expand_archives(root: Path, max_loops: int = 2):
    """Recursively expand nested tar/tgz/tar.gz files within root up to max_loops times."""
    exts = ('.tar', '.tgz', '.tar.gz')
    for _ in range(max_loops):
        tars = [p for p in root.rglob('*') if p.suffix in exts or p.name.endswith('.tar.gz')]
        new_any = False
        for t in tars:
            try:
                with tarfile.open(t.as_posix(), 'r:*') as tf:
                    outdir = t.parent
                    tf.extractall(outdir.as_posix())
                    new_any = True
            except Exception:
                continue
        if not new_any:
            break


def qc_and_basic(adata: ad.AnnData, mito_prefix='MT-'):
    adata.var['mt'] = adata.var_names.str.upper().str.startswith(mito_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    # basic filters (conservative; adjust as needed)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs['pct_counts_mt'] < 20].copy()
    # respect pre-normalized data
    if not adata.uns.get('log_normalized', False):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    n_top = 3000 if adata.n_vars > 3500 else max(1000, int(adata.n_vars*0.8))
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor='seurat_v3', subset=True)
    except ImportError:
        # fallback if scikit-misc not available
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor='seurat', subset=True)
    except Exception:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor='cell_ranger', subset=True)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)
    try:
        # Harmony integration by sample if available
        import harmonypy
        sc.external.pp.harmony_integrate(adata, key='sample')
        sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=20)
    except Exception:
        pass
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)
    return adata


def find_csv_matrices(root: Path):
    """Find pairs of LogNormalized_Data.csv.gz and optional Cell_Types_Annotations.csv.gz per sample."""
    data = {}
    for p in root.rglob('*_LogNormalized_Data.csv.gz'):
        stem = p.name.replace('_LogNormalized_Data.csv.gz','')
        ann = p.parent / f"{stem}_Cell_Types_Annotations.csv.gz"
        data[p] = ann if ann.exists() else None
    return data


def read_csv_matrix(data_csv: Path, ann_csv: Path|None):
    df = pd.read_csv(data_csv, index_col=0)
    # assume rows are genes, columns are cells
    df.index = df.index.astype(str).str.upper()
    df = df[~df.index.duplicated(keep='first')]
    # use sparse float32 to save memory; transpose to cells x genes
    X = sp.csr_matrix(df.values.T.astype(np.float32))
    adata = ad.AnnData(X, var=pd.DataFrame(index=df.index), obs=pd.DataFrame(index=df.columns))
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obs['sample'] = data_csv.stem.split('_')[0]
    adata.uns['log_normalized'] = True
    if ann_csv is not None and ann_csv.exists():
        try:
            ann = pd.read_csv(ann_csv)
            # try merge by barcode/cell id if available
            key = None
            for c in ['barcode','cell','Cell','Cell_ID','cell_id','CellID']:
                if c in ann.columns:
                    key = c; break
            if key is not None:
                ann[key] = ann[key].astype(str)
                obs = adata.obs.copy()
                obs['barcode'] = obs.index.astype(str)
                merged = obs.merge(ann, left_on='barcode', right_on=key, how='left')
                adata.obs = merged.set_index(adata.obs.index)
        except Exception as e:
            print(f"[warn] could not merge annotations for {data_csv.name}: {e}")
    return adata


def main():
    ap = argparse.ArgumentParser(description='Preprocess GSE186344 scRNA: extract/read 10X, QC, integrate, cluster.')
    ap.add_argument('--raw', required=True, help='Path to GSE186344_RAW.tar or extracted folder')
    ap.add_argument('--out', required=True, help='Output h5ad path, e.g., data/processed/scRNA/gse186344.h5ad')
    args = ap.parse_args()

    raw_path = Path(args.raw)
    workdir = None
    if raw_path.is_file() and raw_path.suffix == '.tar':
        tmp = tempfile.mkdtemp(prefix='gse186344_')
        with tarfile.open(raw_path, 'r') as tf:
            tf.extractall(tmp)
        workdir = Path(tmp)
    else:
        workdir = raw_path if raw_path.is_dir() else raw_path.parent

    # expand nested archives if present
    expand_archives(workdir)

    tenx_dirs = find_10x_dirs(workdir)
    tenx_h5 = find_10x_h5(workdir)
    csv_mats = find_csv_matrices(workdir)
    if not tenx_dirs and not tenx_h5 and not csv_mats:
        print('Directory tree sample under', workdir)
        for i, p in enumerate(sorted(workdir.rglob('*'))):
            if i>200: break
            print(p)
        raise SystemExit('No 10X directories, h5, or LogNormalized_Data.csv.gz files found')

    adatas = []
    for d in tenx_dirs:
        try:
            a = read_10x_dir(d)
            a.obs['source_dir'] = d.as_posix()
            adatas.append(a)
        except Exception as e:
            print(f'[warn] skip {d}: {e}')
    for h5 in tenx_h5:
        try:
            a = read_10x_h5(h5)
            a.obs['source_dir'] = h5.as_posix()
            adatas.append(a)
        except Exception as e:
            print(f'[warn] skip {h5}: {e}')
    for data_csv, ann_csv in csv_mats.items():
        try:
            a = read_csv_matrix(data_csv, ann_csv)
            a.obs['source_dir'] = data_csv.as_posix()
            adatas.append(a)
        except Exception as e:
            print(f'[warn] skip {data_csv}: {e}')
    if not adatas:
        raise SystemExit('No readable 10X samples')

    # Concatenate with safer key extraction (avoid positional indexing deprecation), still may be memory-heavy.
    keys = []
    for x in adatas:
        keys.append(x.obs['sample'].iloc[0] if 'sample' in x.obs else 'sample')
    adata = ad.concat(adatas, join='outer', label='batch', keys=keys, fill_value=0)
    adata = qc_and_basic(adata)

    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(outp.as_posix(), compression='gzip')
    print(f'Wrote: {outp}')


if __name__ == '__main__':
    main()
