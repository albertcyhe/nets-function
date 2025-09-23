#!/usr/bin/env python3
"""Build MET500 site-triage prototype using Serpin/ΔFM/THBS1/IL6-STAT3/CoOption features.

Outputs
-------
- results/tables/site_triage_metrics.tsv         (AUC, Brier, etc.)
- results/tables/site_triage_coefficients.tsv    (final model coefficients)
- results/tables/site_triage_predictions.tsv     (out-of-fold probabilities)
- results/tables/site_triage_dca.tsv             (decision-curve values)
- results/figures/A4_site_triage_roc_calib_dca.pdf (combined visual)
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import gzip
import json

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler


SERPIN_CORE = ["SERPINA1", "SERPINA3", "SERPINI1"]
SERPIN_EXTENDED = SERPIN_CORE + ["SERPINB1", "SERPINB2", "SERPINB5", "SERPINB6", "SERPINB8", "SLPI", "PI3", "A2M"]
NET_F_GENES = ["ELANE", "PRTN3", "CTSG"]
IL6_STAT3 = ["IL6", "STAT3", "SOCS3", "JAK1", "JAK2", "IL6ST", "CEBPB", "FOS"]
COOPTION = ["GFAP", "ITGA6", "ITGB1", "LAMC1", "COL4A1", "VEGFA", "ANGPT1", "ANGPT2", "S100B"]
NULL_INFLAMMATORY = ["RPS3", "RPS6", "RPLP0", "RPL11", "EEF1A1", "GAPDH", "ACTB"]
NULL_COOPTION = ["ALDOA", "ENO1", "PGK1", "GPI", "LDHA", "PKM", "TPI1"]


@dataclass
class Inputs:
    expr_path: Path
    clin_path: Path
    gtf_path: Path
    net_m_path: Path
    out_dir: Path


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--expr",
        type=Path,
        default=Path("data/raw/transcriptomics/MET500/M.mx.log2.txt.gz"),
    )
    ap.add_argument(
        "--clin",
        type=Path,
        default=Path("data/raw/transcriptomics/MET500/M.meta.plus.txt"),
    )
    ap.add_argument(
        "--gtf",
        type=Path,
        default=Path("resources/annotation/gencode/gencode.v43.basic.annotation.gtf.gz"),
    )
    ap.add_argument(
        "--net-m",
        type=Path,
        default=Path("resources/modules/net_m_v2.tsv"),
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results"),
    )
    ap.add_argument(
        "--log-dir",
        type=Path,
        default=Path("logs/transcriptomics"),
    )
    ap.add_argument(
        "--variants",
        type=str,
        default="baseline,serpin_extended,null_inflammation",
        help="Comma-separated variant names to evaluate (default: baseline,serpin_extended,null_inflammation)",
    )
    ap.add_argument(
        "--folds", type=int, default=5,
        help="Number of stratified CV folds (default: 5)",
    )
    ap.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed for CV splits.",
    )
    return ap.parse_args()


def read_gtf_mapping(gtf_path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = fields[8]
            info: Dict[str, str] = {}
            for part in attributes.strip().split(";"):
                part = part.strip()
                if part:
                    key, value = part.split(" ", 1)
                    info[key] = value.strip('"')
            gene_id = info.get("gene_id", "").split(".")[0]
            gene_name = info.get("gene_name")
            if gene_id and gene_name:
                mapping[gene_id] = gene_name
    return mapping


def load_expression(expr_path: Path, mapping: Dict[str, str]) -> pd.DataFrame:
    df = pd.read_csv(expr_path, sep="\t", index_col=0)
    gene_ids = df.index.to_series().str.split(".").str[0]
    symbols = gene_ids.map(mapping)
    df.insert(0, "gene_symbol", symbols)
    df = df.dropna(subset=["gene_symbol"])
    df = df.set_index("gene_symbol")
    expr = df.groupby(level=0).mean()
    expr.index = expr.index.str.upper()
    return expr


def module_score(expr: pd.DataFrame, genes: Iterable[str]) -> pd.Series:
    genes = [g for g in genes if g in expr.index]
    if not genes:
        return pd.Series(np.nan, index=expr.columns)
    return expr.loc[genes].mean(axis=0)


def zscore(series: pd.Series) -> pd.Series:
    vals = series.astype(float)
    std = vals.std(ddof=0)
    if std == 0 or np.isnan(std):
        return pd.Series(np.zeros_like(vals), index=series.index)
    return (vals - vals.mean()) / std


def decision_curve(y_true: np.ndarray, probs: np.ndarray, thresholds: np.ndarray) -> pd.DataFrame:
    n = len(y_true)
    prevalence = y_true.mean()
    rows = []
    for thr in thresholds:
        preds = probs >= thr
        tp = np.logical_and(preds, y_true == 1).sum()
        fp = np.logical_and(preds, y_true == 0).sum()
        nb_model = (tp / n) - (fp / n) * (thr / (1 - thr))
        nb_all = prevalence - (1 - prevalence) * (thr / (1 - thr))
        rows.append({
            "threshold": thr,
            "net_benefit_model": nb_model,
            "net_benefit_all": nb_all,
            "net_benefit_none": 0.0,
        })
    return pd.DataFrame(rows)


def load_net_m_genes(path: Path) -> List[str]:
    genes = pd.read_csv(path, header=None)[0].astype(str).tolist()
    return [g.upper() for g in genes if isinstance(g, str)]


def build_features(
    expr: pd.DataFrame,
    serpin_genes: Sequence[str],
    net_f_genes: Sequence[str],
    net_m: Sequence[str],
    il6_genes: Sequence[str],
    cooption_genes: Sequence[str],
) -> pd.DataFrame:
    serpin = module_score(expr, [g.upper() for g in serpin_genes])
    net_f = module_score(expr, [g.upper() for g in net_f_genes])
    net_m_score = module_score(expr, [g.upper() for g in net_m])
    il6 = module_score(expr, [g.upper() for g in il6_genes])
    coop = module_score(expr, [g.upper() for g in cooption_genes])
    thbs1 = expr.loc["THBS1"] if "THBS1" in expr.index else pd.Series(np.nan, index=expr.columns, name="THBS1")

    df = pd.DataFrame({
        "Serpin_score": serpin,
        "NetF_score": net_f,
        "NetM_score": net_m_score,
        "THBS1_expression": thbs1,
        "IL6_STAT3_score": il6,
        "CoOption_score": coop,
    })
    df["DeltaFM_RNA"] = zscore(df["NetF_score"]) - zscore(df["NetM_score"])
    df["Serpin_z"] = zscore(df["Serpin_score"])
    df["IL6_STAT3_z"] = zscore(df["IL6_STAT3_score"])
    df["CoOption_z"] = zscore(df["CoOption_score"])
    df["DeltaFM_z_neg"] = -zscore(df["DeltaFM_RNA"])
    df["THBS1_z_neg"] = -zscore(df["THBS1_expression"])
    features = df[["Serpin_z", "DeltaFM_z_neg", "THBS1_z_neg", "IL6_STAT3_z", "CoOption_z"]]
    return features


def fit_site_model(X: pd.DataFrame, y: pd.Series, folds: int, seed: int) -> Dict[str, object]:
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    preds = np.zeros_like(y, dtype=float)
    coefs = []
    intercepts = []
    for train_idx, test_idx in skf.split(X, y):
        scaler = StandardScaler()
        clf = LogisticRegression(max_iter=1000)
        X_train = scaler.fit_transform(X.iloc[train_idx])
        X_test = scaler.transform(X.iloc[test_idx])
        clf.fit(X_train, y.iloc[train_idx])
        preds[test_idx] = clf.predict_proba(X_test)[:, 1]
        coefs.append(clf.coef_[0])
        intercepts.append(clf.intercept_[0])

    final_scaler = StandardScaler().fit(X)
    final_model = LogisticRegression(max_iter=1000)
    final_model.fit(final_scaler.transform(X), y)

    return {
        "preds": preds,
        "coefs": np.array(coefs),
        "intercepts": np.array(intercepts),
        "final_model": final_model,
        "final_scaler": final_scaler,
    }


def compute_metrics(y: pd.Series, preds: np.ndarray) -> Dict[str, float]:
    metrics = {
        "roc_auc": float(roc_auc_score(y, preds)),
        "average_precision": float(average_precision_score(y, preds)),
        "brier_score": float(brier_score_loss(y, preds)),
    }
    # Balanced accuracy at 0.5 threshold
    labels = (preds >= 0.5).astype(int)
    tp = np.sum((labels == 1) & (y == 1))
    tn = np.sum((labels == 0) & (y == 0))
    fp = np.sum((labels == 1) & (y == 0))
    fn = np.sum((labels == 0) & (y == 1))
    sens = tp / (tp + fn) if (tp + fn) else np.nan
    spec = tn / (tn + fp) if (tn + fp) else np.nan
    metrics.update({
        "sensitivity@0.5": float(sens),
        "specificity@0.5": float(spec),
    })
    return metrics


def write_table(path: Path, data: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(path, sep="\t", index=False)


def plot_panels(out_path: Path, y: pd.Series, preds: np.ndarray, dca: pd.DataFrame) -> None:
    fpr, tpr, _ = roc_curve(y, preds)
    prob_true, prob_pred = calibration_curve(y, preds, n_bins=10, strategy="quantile")
    precision, recall, _ = precision_recall_curve(y, preds)

    sns.set_context("talk")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # ROC
    ax = axes[0]
    ax.plot(fpr, tpr, label=f"AUC = {roc_auc_score(y, preds):.2f}")
    ax.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=1)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC Curve")
    ax.legend(loc="lower right")

    # Calibration
    ax = axes[1]
    ax.plot(prob_pred, prob_true, marker="o", label="Model")
    ax.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=1)
    ax.set_xlabel("Predicted risk")
    ax.set_ylabel("Observed risk")
    ax.set_title("Calibration (quantile bins)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend()

    # Decision curve
    ax = axes[2]
    ax.plot(dca["threshold"], dca["net_benefit_model"], label="Model", color="C0")
    all_key = "net_benefit_all" if "net_benefit_all" in dca.columns else "net_benefit_treat_all"
    none_key = "net_benefit_none" if "net_benefit_none" in dca.columns else "net_benefit_treat_none"
    ax.plot(dca["threshold"], dca[all_key], label="Treat-all", color="C1", linestyle="--")
    ax.plot(dca["threshold"], dca[none_key], label="Treat-none", color="black", linestyle=":")
    ax.set_xlabel("Threshold probability")
    ax.set_ylabel("Net benefit")
    ax.set_title("Decision Curve Analysis")
    ax.legend()

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def write_log(log_path: Path, metrics: Dict[str, float], thresholds: Iterable[float]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "metrics": metrics,
        "thresholds": list(thresholds),
    }
    log_path.write_text(json.dumps(payload, indent=2))


def main() -> None:
    args = parse_args()
    inputs = Inputs(
        expr_path=args.expr,
        clin_path=args.clin,
        gtf_path=args.gtf,
        net_m_path=args.net_m,
        out_dir=args.out_dir,
    )

    mapping = read_gtf_mapping(inputs.gtf_path)
    expr = load_expression(inputs.expr_path, mapping)
    clin = pd.read_csv(inputs.clin_path, sep="\t")
    clin = clin.rename(columns={clin.columns[0]: "Sample_id"})
    clin = clin.set_index("Sample_id")

    common = expr.columns.intersection(clin.index)
    expr = expr.loc[:, common]
    clin = clin.loc[common]

    label = clin["biopsy_tissue"].str.lower().eq("brain").astype(int)

    net_m = load_net_m_genes(inputs.net_m_path)
    variant_names = [v.strip() for v in args.variants.split(",") if v.strip()]
    if not variant_names:
        variant_names = ["baseline"]

    variant_specs: Dict[str, Dict[str, Sequence[str]]] = {
        "baseline": {
            "serpin": SERPIN_CORE,
            "net_f": NET_F_GENES,
            "il6": IL6_STAT3,
            "cooption": COOPTION,
            "label": "Core serpins + IL6/CoOption",
        },
        "serpin_extended": {
            "serpin": SERPIN_EXTENDED,
            "net_f": NET_F_GENES,
            "il6": IL6_STAT3,
            "cooption": COOPTION,
            "label": "Extended serpins (+SERPINB/SLPI/A2M)",
        },
        "null_inflammation": {
            "serpin": SERPIN_CORE,
            "net_f": NET_F_GENES,
            "il6": NULL_INFLAMMATORY,
            "cooption": NULL_COOPTION,
            "label": "Null inflammatory proxies (ribosomal/glycolytic)",
        },
    }

    metrics_rows = []
    coef_rows = []
    pred_rows = []

    baseline_metrics: Dict[str, float] | None = None
    thresholds = np.arange(0.05, 0.91, 0.05)

    out_tables = inputs.out_dir / "tables"
    out_figs = inputs.out_dir / "figures"

    for idx, name in enumerate(variant_names):
        if name not in variant_specs:
            raise ValueError(f"Unknown variant '{name}'. Choices: {sorted(variant_specs)}")
        spec = variant_specs[name]
        feats = build_features(
            expr,
            spec["serpin"],
            spec["net_f"],
            net_m,
            spec["il6"],
            spec["cooption"],
        )
        data = feats.join(label.rename("is_brain"))
        data.index.name = "sample"
        data = data.dropna()

        X = data[["Serpin_z", "DeltaFM_z_neg", "THBS1_z_neg", "IL6_STAT3_z", "CoOption_z"]]
        y = data["is_brain"].astype(int)

        model = fit_site_model(X, y, args.folds, args.random_state)
        preds = model["preds"]
        metrics = compute_metrics(y, preds)
        if baseline_metrics is None:
            baseline_metrics = metrics.copy()

        dca = decision_curve(y.values, preds, thresholds)
        dca_out = dca.rename(columns={
            "net_benefit_model": "net_benefit_model",
            "net_benefit_all": "net_benefit_treat_all",
            "net_benefit_none": "net_benefit_treat_none",
        })

        metrics_rows.extend(
            {
                "variant": name,
                "variant_label": spec["label"],
                "metric": key,
                "value": val,
                "delta_vs_baseline": (val - baseline_metrics[key]) if baseline_metrics else np.nan,
            }
            for key, val in metrics.items()
        )

        coef = pd.DataFrame(model["final_model"].coef_, columns=X.columns)
        coef.insert(0, "intercept", model["final_model"].intercept_)
        coef["variant"] = name
        coef_rows.append(coef)

        preds_df = data.reset_index()[["sample", "is_brain"]]
        preds_df["prob_brain"] = preds
        preds_df["variant"] = name
        pred_rows.append(preds_df)

        if idx == 0:
            # Preserve legacy baseline outputs
            metrics_df = pd.DataFrame([{"metric": k, "value": v} for k, v in metrics.items()])
            write_table(out_tables / "site_triage_metrics.tsv", metrics_df)

            coef_out = coef.copy()
            coef_out.insert(0, "variant", coef_out.pop("variant"))
            write_table(out_tables / "site_triage_coefficients.tsv", coef_out.drop(columns=["variant"]))

            base_preds = preds_df.drop(columns=["variant"])
            write_table(out_tables / "site_triage_predictions.tsv", base_preds)

            write_table(out_tables / "site_triage_dca.tsv", dca_out)
            plot_panels(out_figs / "A4_site_triage_roc_calib_dca.pdf", y, preds, dca_out)

            log_path = args.log_dir / f"met500_site_triage_{args.random_state}.json"
            write_log(log_path, metrics, thresholds)

        write_table(out_tables / f"site_triage_dca_{name}.tsv", dca_out.assign(variant=name))

    # Aggregate outputs
    metrics_table = pd.DataFrame(metrics_rows)
    write_table(out_tables / "site_triage_variant_metrics.tsv", metrics_table)

    coef_table = pd.concat(coef_rows, ignore_index=True)
    write_table(out_tables / "site_triage_variant_coefficients.tsv", coef_table)

    preds_table = pd.concat(pred_rows, ignore_index=True)
    write_table(out_tables / "site_triage_variant_predictions.tsv", preds_table)


if __name__ == "__main__":
    main()
