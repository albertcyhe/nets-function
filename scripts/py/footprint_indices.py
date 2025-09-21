#!/usr/bin/env python3
"""Compute simple NE/PR3 footprint indices from PSM tables.

Heuristic proxy using closed/semi-tryptic searches:
  - Count semi-/non-tryptic PSMs per sample (Number of Enzymatic Termini < 2)
  - Among them, count those with hydrophobic N-terminal residue consistent with
    NE/PR3 (P1 = A/V/I/L/M/F) and non-tryptic context (Prev AA not K/R)
  - Footprint index = log2(motif_count+1) - log2(semi_count+1)

Inputs:
  - data/interim/proteomics_combined/<dataset>/fragger_closed/psm.tsv (preferred)
    or data/interim/proteomics/<dataset>/fragger_closed/psm.tsv
  - data/processed/proteomics/<dataset>/proteo_deltafm.tsv (for joins)

Outputs:
  - results/tables/footprints.tsv – per sample indices
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd


MOTIF_P1 = set(list("AVILMF"))  # NE/PR3-like hydrophobic
CTRL_P1 = set(list("DE"))       # negative-control motif (acidic P1)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", action="append", dest="datasets", required=True)
    ap.add_argument("--combined-root", type=Path, default=Path("data/interim/proteomics_combined"))
    ap.add_argument("--fallback-root", type=Path, default=Path("data/interim/proteomics"))
    ap.add_argument("--processed-root", type=Path, default=Path("data/processed/proteomics"))
    ap.add_argument("--outdir", type=Path, default=Path("results/tables"))
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
        s = s.str.replace("interact-", "", regex=False).str.replace(".pep.xml", "", regex=False)
        return s
    if "Spectrum" in df.columns:
        return df["Spectrum"].astype(str).str.split(".").str[0]
    return pd.Series(["unknown"] * len(df))


def main() -> None:
    args = parse_args()
    out_rows = []

    for ds in args.datasets:
        psm = load_psm(ds, args.combined_root, args.fallback_root)
        psm["Sample"] = derive_sample(psm)

        semicol = None
        for c in psm.columns:
            if c.lower().strip().startswith("number of enzymatic termini"):
                semicol = c
                break
        if semicol is None and "Number of Enzymatic Termini" in psm.columns:
            semicol = "Number of Enzymatic Termini"

        if semicol is not None:
            semi_mask = psm[semicol].fillna(2).astype(float) < 2
        else:
            # fallback: approximate by Prev/Next AA not consistent with trypsin
            semi_mask = (~psm["Prev AA"].astype(str).isin(["K", "R"])) | (~psm["Peptide"].astype(str).str.endswith(("K", "R")))

        nterm = psm["Peptide"].astype(str).str.replace("[^A-Z]", "", regex=True).str[:1]
        prev = psm.get("Prev AA", pd.Series(["?"] * len(psm)))

        motif_mask = semi_mask & nterm.isin(MOTIF_P1) & (~prev.astype(str).isin(["K", "R"]))
        ctrl_mask = semi_mask & nterm.isin(CTRL_P1) & (~prev.astype(str).isin(["K", "R"]))

        grp = psm.groupby("Sample")
        semi_count = grp.apply(lambda g: int(semi_mask.loc[g.index].sum()))
        motif_count = grp.apply(lambda g: int(motif_mask.loc[g.index].sum()))
        ctrl_count  = grp.apply(lambda g: int(ctrl_mask.loc[g.index].sum()))
        idx = pd.DataFrame({
            "dataset": ds,
            "sample": semi_count.index,
            "semi_count": semi_count.values,
            "motif_count": motif_count.values,
            "ctrl_motif_count": ctrl_count.values,
        })
        idx["footprint_index"] = np.log2(idx["motif_count"] + 1.0) - np.log2(idx["semi_count"] + 1.0)
        idx["control_footprint_index"] = np.log2(idx["ctrl_motif_count"] + 1.0) - np.log2(idx["semi_count"] + 1.0)

        # join ΔFM if available
        delta_path = args.processed_root / ds / "proteo_deltafm.tsv"
        if delta_path.exists():
            delta = pd.read_csv(delta_path, sep="\t")[['Sample','Proteo_DeltaFM']]
            idx = idx.merge(delta, left_on="sample", right_on="Sample", how="left").drop(columns=["Sample"])

        out_rows.append(idx)

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    pd.concat(out_rows, ignore_index=True).to_csv(outdir / "footprints.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
