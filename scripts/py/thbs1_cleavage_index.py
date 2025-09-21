#!/usr/bin/env python3
"""Compute THBS1 cleavage index and collagen controls from PSM tables.

Definition (per sample):
  cleave_idx = log2( (# THBS1 semi/non-tryptic PSMs + 1) / (# all THBS1 PSMs + 1) )

Controls: same index for COL1A1 and COL4A1 (non-NE/PR3 preferred substrates)

Inputs:
  - data/interim/proteomics_combined/<dataset>/fragger_closed/psm.tsv (preferred)
  - fallback: data/interim/proteomics/<dataset>/fragger_closed/psm.tsv
  - optional: data/processed/proteomics/<dataset>/proteo_deltafm.tsv (to join ΔFM)
  - results/tables/serpin_scores.tsv (to join Serpin scores)

Output:
  - results/tables/thbs1_cleave_idx.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


GENES = ["THBS1", "COL1A1", "COL4A1"]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", action="append", dest="datasets", required=True)
    ap.add_argument("--combined-root", type=Path, default=Path("data/interim/proteomics_combined"))
    ap.add_argument("--fallback-root", type=Path, default=Path("data/interim/proteomics"))
    ap.add_argument("--processed-root", type=Path, default=Path("data/processed/proteomics"))
    ap.add_argument("--serpin-scores", type=Path, default=Path("results/tables/serpin_scores.tsv"))
    ap.add_argument("--out", type=Path, default=Path("results/tables/thbs1_cleave_idx.tsv"))
    return ap.parse_args()


def load_psm(dataset: str, combined_root: Path, fallback_root: Path) -> pd.DataFrame:
    path = combined_root / dataset / "fragger_closed" / "psm.tsv"
    if not path.exists():
        path = fallback_root / dataset / "fragger_closed" / "psm.tsv"
    df = pd.read_csv(path, sep="\t")
    return df


def derive_sample(df: pd.DataFrame) -> pd.Series:
    if "Spectrum File" in df.columns and df["Spectrum File"].nunique() > 1:
        s = df["Spectrum File"].astype(str)
        return s.str.replace("interact-", "", regex=False).str.replace(".pep.xml", "", regex=False)
    if "Spectrum" in df.columns:
        return df["Spectrum"].astype(str).str.split(".").str[0]
    return pd.Series(["unknown"] * len(df))


def main() -> None:
    args = parse_args()
    scores = pd.read_csv(args.serpin_scores, sep="\t") if args.serpin_scores.exists() else pd.DataFrame()
    rows: List[Dict[str, object]] = []

    for ds in args.datasets:
        psm = load_psm(ds, args.combined_root, args.fallback_root)
        psm["Sample"] = derive_sample(psm)
        semicol = None
        for c in psm.columns:
            if c.lower().strip().startswith("number of enzymatic termini") or c.strip() == "Number of Enzymatic Termini":
                semicol = c
                break
        semi_mask = psm[semicol].fillna(2).astype(float) < 2 if semicol is not None else pd.Series([False]*len(psm))

        # per gene per sample counts
        for gene in GENES:
            gmask = psm.get("Gene", pd.Series(["?"] * len(psm))).astype(str) == gene
            if not gmask.any():
                continue
            df_gene = psm[gmask].copy()
            df_gene = df_gene.assign(is_semi=semi_mask[gmask].astype(int).values)
            tmp = df_gene.groupby("Sample").agg(coverage=("is_semi", "size"), semi=("is_semi", "sum")).reset_index()
            tmp.insert(0, "dataset", ds)
            tmp.rename(columns={"Sample": "sample", "coverage": f"{gene}_coverage", "semi": f"{gene}_semi"}, inplace=True)
            tmp[f"{gene}_cleave_idx"] = np.log2(tmp[f"{gene}_semi"] + 1.0) - np.log2(tmp[f"{gene}_coverage"] + 1.0)
            rows.append(tmp)

    out = None
    if rows:
        out = pd.DataFrame()
        for r in rows:
            out = r if out.empty else out.merge(r, on=["dataset","sample"], how="outer")
        # join ΔFM
        delta_rows = []
        for ds in args.datasets:
            dpath = args.processed_root / ds / "proteo_deltafm.tsv"
            if dpath.exists():
                d = pd.read_csv(dpath, sep="\t")[['Sample','Proteo_DeltaFM']]
                d["dataset"] = ds
                d.rename(columns={"Sample": "sample"}, inplace=True)
                delta_rows.append(d)
        if delta_rows:
            delta = pd.concat(delta_rows, ignore_index=True)
            out = out.merge(delta, on=["dataset","sample"], how="left")
        # join Serpin scores
        if not scores.empty:
            out = out.merge(scores[["dataset","sample","Serpin_score_core","Serpin_score_ext"]], on=["dataset","sample"], how="left")
        # fill zeros for missing counts
        for col in [f"{g}_coverage" for g in GENES] + [f"{g}_semi" for g in GENES]:
            if col in out.columns:
                out[col] = out[col].fillna(0).astype(int)

        args.out.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
