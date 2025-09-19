#!/usr/bin/env python3
import argparse
import gzip
import os
import sys
import pandas as pd


def read_series_matrix(path):
    opener = gzip.open if path.endswith('.gz') else open
    meta = {}
    table_started = False
    rows = []
    with opener(path, 'rt', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('!') and not table_started:
                # collect selected meta lines
                if line.startswith('!Series_platform_id'):
                    meta['Series_platform_id'] = line.split('\t')[-1].strip('"')
                if line.startswith('!Sample_title'):
                    meta['Sample_title'] = [s.strip('"') for s in line.split('\t')[1:]]
                if line.startswith('!Sample_geo_accession'):
                    meta['Sample_geo_accession'] = [s.strip('"') for s in line.split('\t')[1:]]
                if line.startswith('!Sample_characteristics_ch1'):
                    meta.setdefault('Sample_characteristics_ch1', []).append([s.strip('"') for s in line.split('\t')[1:]])
                continue
            if line.startswith('!series_matrix_table_begin'):
                table_started = True
                header = None
                continue
            if line.startswith('!series_matrix_table_end'):
                break
            if table_started:
                parts = line.split('\t')
                if header is None:
                    header = parts
                else:
                    rows.append(parts)
    if not rows:
        # Fallback: re-read capturing table lines verbatim and parse with csv
        import csv
        opener = gzip.open if path.endswith('.gz') else open
        table_lines = []
        with opener(path, 'rt', encoding='utf-8', errors='ignore') as fh2:
            cap = False
            for ln in fh2:
                if ln.startswith('!series_matrix_table_begin'):
                    cap = True
                    continue
                if ln.startswith('!series_matrix_table_end'):
                    break
                if cap:
                    table_lines.append(ln.rstrip('\n'))
        if not table_lines:
            raise RuntimeError('No expression table found in series matrix')
        reader = csv.reader(table_lines, delimiter='\t', quotechar='"')
        rows2 = list(reader)
        header = rows2[0]
        rows = rows2[1:]
    df = pd.DataFrame(rows, columns=header)
    df = df.set_index(df.columns[0])
    # coerce to numeric where possible
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    return df, meta


def read_gpl_mapping(gpl_path):
    opener = gzip.open if gpl_path.endswith('.gz') else open
    header = None
    rows = []
    with opener(gpl_path, 'rt', encoding='utf-8', errors='ignore') as fh:
        table_started = False
        for line in fh:
            line = line.rstrip('\n')
            if not table_started:
                if line.startswith('!platform_table_begin'):
                    table_started = True
                    # next line will be header
                    header = None
                continue
            # after table begin
            if header is None:
                header = line.split('\t')
                continue
            if line.startswith('!platform_table_end'):
                break
            parts = line.split('\t')
            # pad/truncate to header length
            if len(parts) < len(header):
                parts = parts + [''] * (len(header) - len(parts))
            elif len(parts) > len(header):
                parts = parts[:len(header)]
            rows.append(parts)
    gpl = pd.DataFrame(rows, columns=header)
    # heuristics: prefer official symbol columns
    candidate_cols = [
        'Gene Symbol', 'Gene symbol', 'Symbol', 'GENE_SYMBOL', 'Gene ID',
        'Gene Symbol;','Gene Symbol\n','ENTREZ_GENE_ID','ORF','Associated Gene Name'
    ]
    symbol_col = None
    for c in candidate_cols:
        if c in gpl.columns:
            symbol_col = c
            break
    if symbol_col is None:
        # try to parse from 'ID' if pipe-delimited
        symbol_col = 'Gene Symbol' if 'Gene Symbol' in gpl.columns else None
    id_col = 'ID' if 'ID' in gpl.columns else gpl.columns[0]
    gpl = gpl[[id_col] + ([symbol_col] if symbol_col else [])].copy()
    gpl.rename(columns={id_col: 'PROBE', symbol_col: 'SYMBOL' if symbol_col else 'SYMBOL'}, inplace=True)
    if 'SYMBOL' not in gpl.columns:
        gpl['SYMBOL'] = ''
    # clean symbols (split on delimiters, take first)
    gpl['SYMBOL'] = gpl['SYMBOL'].astype(str).str.replace(' /// ', ';').str.replace(' // ', ';').str.replace(' | ', ';')
    gpl['SYMBOL'] = gpl['SYMBOL'].str.split(';').str[0].str.strip()
    gpl = gpl[['PROBE', 'SYMBOL']]
    return gpl


def collapse_probes(expr, mapping, method='max'):
    df = expr.copy()
    df.index.name = 'PROBE'
    df = df.merge(mapping, left_index=True, right_on='PROBE', how='left')
    df['SYMBOL'] = df['SYMBOL'].fillna('')
    df = df[df['SYMBOL'] != '']
    value_cols = [c for c in df.columns if c not in ('PROBE', 'SYMBOL')]
    if method == 'max':
        agg = df.groupby('SYMBOL')[value_cols].max()
    elif method == 'median':
        agg = df.groupby('SYMBOL')[value_cols].median()
    else:
        agg = df.groupby('SYMBOL')[value_cols].mean()
    return agg


def write_pheno(meta, out_pheno):
    pheno = pd.DataFrame({'sample_id': meta.get('Sample_geo_accession', [])})
    if 'Sample_title' in meta:
        pheno['title'] = meta['Sample_title']
    # unpack characteristics
    chars = meta.get('Sample_characteristics_ch1', [])
    for i, arr in enumerate(chars):
        key = f'char{i+1}'
        pheno[key] = arr
    pheno.to_csv(out_pheno, sep='\t', index=False)


def main():
    ap = argparse.ArgumentParser(description='Parse GEO series matrix, map probes via GPL, export expr/pheno TSVs.')
    ap.add_argument('--series', required=True, help='Path to GSE*_series_matrix.txt[.gz]')
    ap.add_argument('--gpl', required=True, help='Path to GPL*.annot[.gz]')
    ap.add_argument('--dataset', required=True, help='Dataset name (e.g., GSE125989)')
    ap.add_argument('--outdir', required=True, help='Output directory under data/processed/<dataset>')
    ap.add_argument('--collapse', default='max', choices=['max','median','mean'])
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    expr_raw, meta = read_series_matrix(args.series)
    mapping = read_gpl_mapping(args.gpl)
    expr = collapse_probes(expr_raw, mapping, method=args.collapse)

    out_expr = os.path.join(args.outdir, f'{args.dataset}.expr.tsv.gz')
    out_pheno = os.path.join(args.outdir, f'{args.dataset}.pheno.tsv')
    out_platform = os.path.join(args.outdir, f'{args.dataset}.platform.tsv')

    expr.to_csv(out_expr, sep='\t', index=True, header=True, compression='gzip')
    write_pheno(meta, out_pheno)
    # minimal platform info
    plat = pd.DataFrame({'Series_platform_id':[meta.get('Series_platform_id','NA')], 'series_file':[os.path.basename(args.series)], 'gpl_file':[os.path.basename(args.gpl)]})
    plat.to_csv(out_platform, sep='\t', index=False)
    print(f'Wrote: {out_expr}\nWrote: {out_pheno}\nWrote: {out_platform}')


if __name__ == '__main__':
    sys.exit(main())
