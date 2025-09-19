#!/usr/bin/env python3
import argparse
import os
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description='Process GSE184869 Excel (log2 TMM CPM protein-coding).')
    ap.add_argument('--excel', required=True, help='Path to GSE184869 Excel file')
    ap.add_argument('--dataset', default='GSE184869')
    ap.add_argument('--outdir', required=True, help='Output directory under data/processed/<dataset>')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    # Read Excel; auto-detect sheet and format
    xls = pd.ExcelFile(args.excel)
    # pick the first sheet by default
    sheet = xls.sheet_names[0]
    df = pd.read_excel(args.excel, sheet_name=sheet)

    # heuristics: gene column detection
    gene_cols = [c for c in df.columns if str(c).lower() in ('gene','genes','symbol','gene_name','hgnc_symbol','ensembl','ensembl_gene_id')]
    if gene_cols:
        gene_col = gene_cols[0]
    else:
        gene_col = df.columns[0]
    df.rename(columns={gene_col: 'SYMBOL'}, inplace=True)
    df['SYMBOL'] = df['SYMBOL'].astype(str)
    # drop empty gene rows
    df = df[df['SYMBOL'].notna() & (df['SYMBOL']!='')]
    # set index and keep numeric columns only for expression
    df = df.set_index('SYMBOL')
    # Coerce expression to numeric
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    # drop all-NA rows
    df = df.dropna(how='all')

    out_expr = os.path.join(args.outdir, f'{args.dataset}.expr.tsv.gz')
    out_pheno = os.path.join(args.outdir, f'{args.dataset}.pheno.tsv')
    out_platform = os.path.join(args.outdir, f'{args.dataset}.platform.tsv')

    df.to_csv(out_expr, sep='\t', index=True, header=True, compression='gzip')

    # Minimal pheno from column names
    pheno = pd.DataFrame({'sample_id': list(df.columns)})
    pheno.to_csv(out_pheno, sep='\t', index=False)

    plat = pd.DataFrame({'platform':['RNA-seq'], 'normalization':['log2 TMM CPM (provided)']})
    plat.to_csv(out_platform, sep='\t', index=False)

    print(f'Wrote: {out_expr}\nWrote: {out_pheno}\nWrote: {out_platform}')


if __name__ == '__main__':
    main()

