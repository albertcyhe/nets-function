#!/usr/bin/env python3
"""Convert GEO bigWig tracks into gene-level coverage matrices.

Designed for RNA-seq GEO series that distribute per-sample stranded bigWig
tracks (e.g. GSE96860). Coverage is summarised over gene bodies defined by a
reference GTF, aggregating plus/minus strands prior to export.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import pyBigWig


logger = logging.getLogger(__name__)


@dataclass
class GeneWindow:
    chrom: str
    start: int  # 0-based inclusive
    end: int    # 0-based exclusive
    gene_id: str
    gene_name: str


def parse_attributes(attr: str) -> Dict[str, str]:
    fields = [x.strip() for x in attr.split(';') if x.strip()]
    out: Dict[str, str] = {}
    for field in fields:
        if ' ' not in field:
            continue
        key, value = field.split(' ', 1)
        out[key] = value.strip('"')
    return out


def load_genes(gtf_path: Path, chromosomes: Iterable[str] | None = None) -> List[GeneWindow]:
    chrom_filter = set(chromosomes) if chromosomes else None
    opener = gzip.open if gtf_path.suffix == '.gz' else open
    genes: List[GeneWindow] = []
    with opener(gtf_path, 'rt', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, _, feature, start, end, _, _, _, attr = parts
            if feature != 'gene':
                continue
            if chrom_filter and chrom not in chrom_filter:
                continue
            start0 = max(int(start) - 1, 0)
            end0 = int(end)
            attrs = parse_attributes(attr)
            gene_id = attrs.get('gene_id', '')
            gene_name = attrs.get('gene_name', gene_id)
            if not gene_id and not gene_name:
                continue
            genes.append(GeneWindow(chrom, start0, end0, gene_id, gene_name))
    if not genes:
        raise RuntimeError(f'No gene entries parsed from {gtf_path}')
    logger.info('Loaded %d genes from %s', len(genes), gtf_path)
    return genes


def summarise_bigwig(pair: Tuple[Path, Path], genes: List[GeneWindow]) -> pd.Series:
    plus_path, minus_path = pair
    with pyBigWig.open(str(plus_path)) as bw_plus, pyBigWig.open(str(minus_path)) as bw_minus:
        chrom_sizes = bw_plus.chroms()
        values = []
        for gene in genes:
            chrom_len = chrom_sizes.get(gene.chrom)
            if chrom_len is None:
                values.append(np.nan)
                continue
            start = min(gene.start, chrom_len - 1)
            end = min(gene.end, chrom_len)
            if start >= end:
                values.append(np.nan)
                continue
            plus = bw_plus.stats(gene.chrom, start, end, type='sum')[0]
            minus = bw_minus.stats(gene.chrom, start, end, type='sum')[0]
            v = 0.0
            if plus is not None and not np.isnan(plus):
                v += plus
            if minus is not None and not np.isnan(minus):
                v += minus
            values.append(v)
    data = pd.Series(values, index=[g.gene_name or g.gene_id for g in genes], dtype=float)
    return data


def collect_bigwig_pairs(directory: Path) -> Dict[str, Tuple[Path, Path]]:
    plus_files = sorted(directory.glob('*_PlusStrand.bw'))
    pattern = re.compile(r'GSE\d+_(.+?)_RNAseq_(.+)_PlusStrand', re.IGNORECASE)
    pairs: Dict[str, Tuple[Path, Path]] = {}
    for plus in plus_files:
        match = pattern.search(plus.stem)
        if not match:
            logger.debug('Skipping unrecognised bigwig: %s', plus.name)
            continue
        cell, condition = match.groups()
        sample_id = f'{cell}_{condition}'.replace('-', '').replace(' ', '')
        minus = plus.with_name(plus.name.replace('PlusStrand', 'MinusStrand'))
        if not minus.exists():
            logger.warning('Minus strand file missing for %s', plus.name)
            continue
        pairs[sample_id] = (plus, minus)
    if not pairs:
        raise RuntimeError(f'No plus/minus bigwig pairs detected under {directory}')
    logger.info('Detected %d bigwig pairs', len(pairs))
    return pairs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bigwig-dir', required=True, type=Path, help='Directory containing *_PlusStrand.bw and *_MinusStrand.bw files')
    ap.add_argument('--gtf', required=True, type=Path, help='Reference gene annotation (GTF, possibly gzipped)')
    ap.add_argument('--dataset', required=True, help='Dataset identifier (e.g. GSE96860)')
    ap.add_argument('--output-dir', required=True, type=Path, help='Destination inside data/processed/<dataset>/')
    ap.add_argument('--chromosomes', nargs='*', help='Optional list of chromosomes to keep (e.g. chr1 chr2 ... chrX)')
    ap.add_argument('--log-level', default='INFO')
    args = ap.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format='%(levelname)s: %(message)s')

    bigwig_dir = args.bigwig_dir
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    genes = load_genes(args.gtf, args.chromosomes)
    pairs = collect_bigwig_pairs(bigwig_dir)

    matrix = {}
    pheno_records = []
    for sample_id, pair in pairs.items():
        logger.info('Summarising %s', sample_id)
        expr = summarise_bigwig(pair, genes)
        matrix[sample_id] = expr
        cell, condition = sample_id.split('_', 1)
        pheno_records.append({'sample_id': sample_id, 'cell_line': cell, 'condition': condition})

    expr_df = pd.DataFrame(matrix).fillna(0.0)
    expr_df.index.name = 'gene'
    expr_path = output_dir / f'{args.dataset}.expr.tsv.gz'
    expr_df.to_csv(expr_path, sep='\t', compression='gzip')

    pheno_df = pd.DataFrame(pheno_records)
    pheno_path = output_dir / f'{args.dataset}.pheno.tsv'
    pheno_df.to_csv(pheno_path, sep='\t', index=False)

    platform_path = output_dir / f'{args.dataset}.platform.tsv'
    platform_df = pd.DataFrame([
        {'Series_platform_id': 'RNAseq_bigwig', 'gtf': args.gtf.name, 'notes': 'Summarised from strand-specific bigWig coverage'}
    ])
    platform_df.to_csv(platform_path, sep='\t', index=False)

    logger.info('Expression matrix written to %s', expr_path)


if __name__ == '__main__':
    main()
