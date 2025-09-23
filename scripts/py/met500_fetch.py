#!/usr/bin/env python3
"""Download MET500 RNA-seq expression and clinical metadata from UCSC Xena."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Tuple

import pandas as pd
import requests

XENA_HOST = "https://ucscpublic.xenahubs.net"
EXPRESSION_DATASET = "MET500/geneExpression/M.mx.log2.txt.gz"  # log2 TPM matrix, genes x samples
CLINICAL_DATASET = "MET500/geneExpression/M.meta.plus.txt"  # clinical/phenotype annotations


def download_dataset(dataset: str, out_path: Path, host: str = XENA_HOST) -> None:
    url = f"{host}/download/{dataset}"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with out_path.open("wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)


def load_dataset(path: Path) -> pd.DataFrame:
    compression = "gzip" if path.suffix == ".gz" else None
    return pd.read_table(path, sep="\t", index_col=0, compression=compression)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-dir", type=Path, default=Path("data/raw/transcriptomics/MET500"))
    ap.add_argument("--processed-dir", type=Path, default=Path("data/processed/transcriptomics/MET500"))
    ap.add_argument(
        "--exclude-tissues",
        type=lambda s: [x.strip() for x in s.split(",") if x.strip()],
        default=["brain", "csf"],
        help="Comma-separated list of tissues (case-insensitive) to exclude.",
    )
    ap.add_argument(
        "--host",
        default=XENA_HOST,
        help="Xena host providing MET500 datasets (default: ucscpublic hub)",
    )
    ap.add_argument(
        "--expr-dataset",
        default=EXPRESSION_DATASET,
        help="Expression dataset path relative to /download/",
    )
    ap.add_argument(
        "--clinical-dataset",
        default=CLINICAL_DATASET,
        help="Clinical dataset path relative to /download/",
    )
    args = ap.parse_args()

    raw_dir = args.raw_dir
    processed_dir = args.processed_dir

    expr_raw_path = raw_dir / Path(args.expr_dataset).name
    clin_raw_path = raw_dir / Path(args.clinical_dataset).name

    print(f"Downloading expression to {expr_raw_path} ...", file=sys.stderr)
    download_dataset(args.expr_dataset, expr_raw_path, host=args.host)

    print(f"Downloading clinical to {clin_raw_path} ...", file=sys.stderr)
    download_dataset(args.clinical_dataset, clin_raw_path, host=args.host)

    print("Parsing expression ...", file=sys.stderr)
    expr = load_dataset(expr_raw_path)

    print("Parsing clinical ...", file=sys.stderr)
    clin = load_dataset(clin_raw_path)

    exclude_terms = {x.lower() for x in args.exclude_tissues}
    biopsy_lower = clin.get("biopsy_tissue", pd.Series(index=clin.index, dtype="object")).astype(str).str.lower()
    tissue_lower = clin.get("tissue", pd.Series(index=clin.index, dtype="object")).astype(str).str.lower()
    exclude_mask = biopsy_lower.isin(exclude_terms) | tissue_lower.isin(exclude_terms)

    kept = clin.loc[~exclude_mask].copy()
    dropped = clin.loc[exclude_mask].copy()

    if dropped.empty:
        print("No samples matched the exclusion criteria.", file=sys.stderr)
    else:
        print(
            f"Excluded {len(dropped)} samples ({', '.join(sorted(set(dropped['biopsy_tissue'].dropna().astype(str))))})",
            file=sys.stderr,
        )

    shared = [s for s in kept.index if s in expr.columns]
    missing_in_expr = sorted(set(kept.index) - set(shared))
    if missing_in_expr:
        print(
            f"Warning: {len(missing_in_expr)} samples present in clinical but missing from expression.",
            file=sys.stderr,
        )
    expr_filtered = expr.loc[:, shared].copy()
    kept = kept.loc[shared]

    processed_dir.mkdir(parents=True, exist_ok=True)

    expr_out = processed_dir / "MET500_expression_log2.nonbrain.tsv.gz"
    clin_out = processed_dir / "MET500_clinical.nonbrain.tsv"
    expr_filtered.to_csv(expr_out, sep="\t", compression="gzip")
    kept.to_csv(clin_out, sep="\t")

    counts_out = processed_dir / "MET500_biopsy_tissue_counts.tsv"
    if "biopsy_tissue" in kept.columns:
        kept["biopsy_tissue"].value_counts().sort_index().to_csv(counts_out, sep="\t", header=["count"])
    else:
        pd.Series(dtype="int64").to_csv(counts_out, sep="\t", header=["count"])

    print(
        f"Saved {expr_filtered.shape[1]} samples across {expr_filtered.shape[0]} genes to {expr_out}",
        file=sys.stderr,
    )
    print(
        f"Clinical annotations (non-brain/CSF) written to {clin_out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
