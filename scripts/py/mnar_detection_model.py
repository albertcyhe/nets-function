#!/usr/bin/env python3
"""MNAR detection models: does Serpin burden predict missingness of F targets?

For each dataset, constructs per-sample detection flags for F (ELANE/PRTN3/CTSG)
and M (NET-M markers) from gene_psm_matrix (log2(count+1)). Fits logistic
regressions:

  logit Pr(detect_F) ~ Serpin_score + log_total_PSMs
  logit Pr(detect_M) ~ Serpin_score + log_total_PSMs   (control)

Outputs:
  - results/tables/mnar_detection_table.tsv – sample-level predictors/responses
  - results/tables/mnar_logit_results.tsv – per-dataset coefficients
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import statsmodels.api as sm


F_TARGETS = ["ELANE", "PRTN3", "CTSG"]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", action="append", dest="datasets", required=True)
    ap.add_argument("--processed-root", type=Path, default=Path("data/processed/proteomics"))
    ap.add_argument("--netm-file", type=Path, default=Path("resources/modules/net_m_proteome.tsv"))
    ap.add_argument("--outdir", type=Path, default=Path("results/tables"))
    return ap.parse_args()


def load_gene_psm(ds_dir: Path) -> pd.DataFrame:
    path = ds_dir / "gene_psm_matrix.tsv"
    mat = pd.read_csv(path, sep="\t", index_col=0)
    counts = np.power(2.0, mat) - 1.0
    return counts


def main() -> None:
    args = parse_args()
    netm = pd.read_csv(args.netm_file, sep="\t")["Gene"].dropna().astype(str).tolist()
    out_rows = []
    result_rows = []

    for ds in args.datasets:
        ds_dir = args.processed_root / ds
        counts = load_gene_psm(ds_dir)
        # Orient gene x sample if necessary
        probe_genes = set(F_TARGETS + netm)
        if len(probe_genes.intersection(counts.columns)) > len(probe_genes.intersection(counts.index)):
            counts = counts.T
        samples = counts.columns

        # predictors
        inhibitors = [g for g in ["SERPINA1", "SERPINA3", "SERPINB1", "SERPINB6", "SERPINB8", "SLPI", "PI3", "ELAFIN", "A2M"] if g in counts.index]
        serpin_score = counts.loc[inhibitors].replace(0, np.nan).apply(np.log2).mean(axis=0) if inhibitors else pd.Series(np.nan, index=samples)
        total_psms = counts.sum(axis=0)

        # responses
        f_present = [g for g in F_TARGETS if g in counts.index]
        m_present = [g for g in netm if g in counts.index]
        detect_f = (counts.loc[f_present] > 0).any(axis=0) if f_present else pd.Series(False, index=samples)
        detect_m = (counts.loc[m_present] > 0).any(axis=0) if m_present else pd.Series(False, index=samples)

        df = pd.DataFrame({
            "dataset": ds,
            "sample": samples,
            "detect_F": detect_f.astype(int).values,
            "detect_M": detect_m.astype(int).values,
            "Serpin_score": serpin_score.values,
            "log_total_PSMs": np.log1p(total_psms.values),
        })
        out_rows.append(df)

        # logistic regressions (per dataset)
        for resp in ("detect_F", "detect_M"):
            # only run if both classes present
            if df[resp].nunique() < 2:
                result_rows.append({"dataset": ds, "response": resp, "n": int(len(df)), "beta_serpin": np.nan, "se": np.nan, "p_value": np.nan})
                continue
            X = sm.add_constant(df[["Serpin_score", "log_total_PSMs"]].fillna(0))
            y = df[resp]
            try:
                model = sm.Logit(y, X).fit(disp=False)
                result_rows.append({
                    "dataset": ds,
                    "response": resp,
                    "n": int(len(df)),
                    "beta_serpin": float(model.params.get("Serpin_score", np.nan)),
                    "se": float(model.bse.get("Serpin_score", np.nan)),
                    "p_value": float(model.pvalues.get("Serpin_score", np.nan)),
                })
            except Exception:
                result_rows.append({"dataset": ds, "response": resp, "n": int(len(df)), "beta_serpin": np.nan, "se": np.nan, "p_value": np.nan})

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    if out_rows:
        pd.concat(out_rows, ignore_index=True).to_csv(outdir / "mnar_detection_table.tsv", sep="\t", index=False)
    pd.DataFrame(result_rows).to_csv(outdir / "mnar_logit_results.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
