#!/usr/bin/env python3
"""Run stratified, robust associations to consolidate the evidence chain.

Models (per stratum: Brain vs NonBrain; and Combined):
  1) Footprint_index ~ Serpin_score_core + log_total_PSMs
  2) THBS1_log_count ~ Serpin_score_core + log_total_PSMs
  3) Proteo_DeltaFM ~ Serpin_score_core + log_total_PSMs

For each, report:
  - Spearman rho/p
  - OLS beta(Serpin_score_core)/p (HC3 robust SE)
  - Robust RLM (Huber) beta and z/p (approx.)
  - BH-adjusted q across the three endpoints within each stratum

Inputs:
  - results/tables/serpin_scores.tsv
  - results/tables/footprints.tsv
  - data/processed/proteomics/<dataset>/proteo_deltafm.tsv
  - data/processed/proteomics/<dataset>/gene_psm_matrix.tsv (for total PSMs)

Output:
  - results/tables/assoc_summary.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--datasets", nargs="+", default=["PXD046330", "PXD005719", "PXD051579"])
    ap.add_argument("--processed-root", type=Path, default=Path("data/processed/proteomics"))
    ap.add_argument("--scores", type=Path, default=Path("results/tables/serpin_scores.tsv"))
    ap.add_argument("--footprints", type=Path, default=Path("results/tables/footprints.tsv"))
    ap.add_argument("--thbs1cleave", type=Path, default=Path("results/tables/thbs1_cleave_idx.tsv"))
    ap.add_argument("--out", type=Path, default=Path("results/tables/assoc_summary.tsv"))
    return ap.parse_args()


def derive_organ(dataset: str, sample: str) -> str:
    s = sample.upper()
    if dataset == "PXD005719":
        return "Brain" if "BR" in s else "NonBrain"
    if dataset in ("PXD046330", "PXD051579"):
        return "NonBrain"
    return "NonBrain"


def load_deltafm(processed_root: Path, datasets: List[str]) -> pd.DataFrame:
    frames = []
    for ds in datasets:
        path = processed_root / ds / "proteo_deltafm.tsv"
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        df = df[["Sample", "Proteo_DeltaFM"] + (["THBS1_log_count"] if "THBS1_log_count" in df.columns else [])]
        df["dataset"] = ds
        frames.append(df.rename(columns={"Sample": "sample"}))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["dataset","sample","Proteo_DeltaFM","THBS1_log_count"])


def load_total_psms(processed_root: Path, datasets: List[str]) -> pd.DataFrame:
    rows = []
    for ds in datasets:
        matp = processed_root / ds / "gene_psm_matrix.tsv"
        if not matp.exists():
            continue
        mat = pd.read_csv(matp, sep="\t", index_col=0)
        counts = np.power(2.0, mat) - 1.0
        # orient genes x samples
        if counts.shape[0] < counts.shape[1] and set(["SERPINA1","SERPINA3","SERPINI1"]).intersection(counts.columns):
            counts = counts.T
        total = counts.sum(axis=0)
        for s, v in total.items():
            rows.append({"dataset": ds, "sample": s, "log_total_PSMs": float(np.log1p(v))})
    return pd.DataFrame(rows)


def bh_adjust(pvals: pd.Series) -> pd.Series:
    p = pvals.copy().astype(float)
    mask = p.notna()
    n = int(mask.sum())
    if n == 0:
        return p
    sub = p[mask]
    order = sub.sort_values().index
    ranks = pd.Series(np.arange(1, n + 1, dtype=float), index=order)
    qsub = (sub.loc[order] * n / ranks).cummin().clip(upper=1.0)
    out = p.copy()
    out.loc[order] = qsub
    return out


def fit_one(y: pd.Series, X: pd.DataFrame) -> Dict[str, float]:
    res: Dict[str, float] = {}
    # Spearman
    v = pd.concat([y, X], axis=1).dropna()
    if len(v) >= 3:
        rho, p = stats.spearmanr(v.iloc[:, 0], v["Serpin_score_core"])  # first col is y
        res.update({"spearman_rho": rho, "spearman_p": p, "n": int(len(v))})
    else:
        res.update({"spearman_rho": np.nan, "spearman_p": np.nan, "n": int(len(v))})
    # OLS HC3
    try:
        ols = sm.OLS(v.iloc[:, 0], sm.add_constant(v[["Serpin_score_core", "log_total_PSMs"]])).fit(cov_type="HC3")
        res.update({
            "ols_beta": float(ols.params.get("Serpin_score_core", np.nan)),
            "ols_se": float(ols.bse.get("Serpin_score_core", np.nan)),
            "ols_p": float(ols.pvalues.get("Serpin_score_core", np.nan)),
            "ols_r2": float(ols.rsquared),
        })
    except Exception:
        res.update({"ols_beta": np.nan, "ols_se": np.nan, "ols_p": np.nan, "ols_r2": np.nan})
    # Robust RLM (Huber)
    try:
        rlm = sm.RLM(v.iloc[:, 0], sm.add_constant(v[["Serpin_score_core", "log_total_PSMs"]]), M=sm.robust.norms.HuberT()).fit()
        beta = float(rlm.params.get("Serpin_score_core", np.nan))
        bse = float(rlm.bse.get("Serpin_score_core", np.nan))
        z = beta / bse if (bse and np.isfinite(bse)) else np.nan
        p = float(2 * (1 - stats.norm.cdf(abs(z)))) if np.isfinite(z) else np.nan
        res.update({"rlm_beta": beta, "rlm_se": bse, "rlm_p": p})
    except Exception:
        res.update({"rlm_beta": np.nan, "rlm_se": np.nan, "rlm_p": np.nan})
    return res


def main() -> None:
    args = parse_args()
    serpin = pd.read_csv(args.scores, sep="\t")
    serpin = serpin[["dataset","sample","Serpin_score_core"]]
    footprints = pd.read_csv(args.footprints, sep="\t")
    delta = load_deltafm(args.processed_root, args.datasets)
    covar = load_total_psms(args.processed_root, args.datasets)
    # THBS1 cleavage index (best-effort unify column)
    cleave = pd.DataFrame()
    if args.thbs1cleave.exists():
        t = pd.read_csv(args.thbs1cleave, sep='\t')
        # prefer unified column; else pick first matching
        if 'THBS1_cleave_idx' not in t.columns:
            cand = [c for c in t.columns if c.startswith('THBS1_cleave_idx')]
            if cand:
                t = t.rename(columns={cand[0]:'THBS1_cleave_idx'})
        cleave = t[['dataset','sample','THBS1_cleave_idx']].copy()

    # Merge master table
    df = serpin.merge(footprints[["dataset","sample","footprint_index"]], on=["dataset","sample"], how="left")
    df = df.merge(delta, on=["dataset","sample"], how="left")
    df = df.merge(covar, on=["dataset","sample"], how="left")
    if not cleave.empty:
        df = df.merge(cleave, on=["dataset","sample"], how="left")
    df["organ"] = [derive_organ(d, s) for d, s in zip(df["dataset"], df["sample"])]

    # Build strata
    strata = ["Combined", "Brain", "NonBrain"]
    rows: List[Dict[str, object]] = []
    for stratum in strata:
        if stratum == "Combined":
            sub = df.copy()
        else:
            sub = df[df["organ"] == stratum].copy()
        if sub.empty:
            continue
        # Endpoints and signs
        endpoints = {
            "footprint_index": "neg",
            "THBS1_log_count": "pos",
            "Proteo_DeltaFM": "neg",
            "THBS1_cleave_idx": "neg",
        }
        for ycol, expect in endpoints.items():
            y = sub[ycol]
            X = sub[["Serpin_score_core", "log_total_PSMs"]]
            res = fit_one(y, X)
            res.update({"stratum": stratum, "endpoint": ycol, "expect": expect})
            rows.append(res)

    out = pd.DataFrame(rows)
    # BH q per stratum
    out["q_spearman"] = np.nan
    out["q_ols"] = np.nan
    out["q_rlm"] = np.nan
    for stratum in out["stratum"].unique():
        mask = out["stratum"] == stratum
        out.loc[mask, "q_spearman"] = bh_adjust(out.loc[mask, "spearman_p"])  # type: ignore
        out.loc[mask, "q_ols"] = bh_adjust(out.loc[mask, "ols_p"])  # type: ignore
        out.loc[mask, "q_rlm"] = bh_adjust(out.loc[mask, "rlm_p"])  # type: ignore

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
