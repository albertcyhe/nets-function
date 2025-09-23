#!/usr/bin/env python3
"""Generate publication-ready figures for liquid site-triage analyses.

Figures
=======
1. ROC comparison for CSF (MSV000089062) and combined plasma (PXD032767 + PXD018301)
2. Regression effect sizes (robust logit + mixed-effects) for the combined plasma model
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import roc_curve, auc

sns.set_context("talk")

RESULTS_ROOT = Path("results/tables")
FIG_ROOT = Path("results/figures")
FIG_ROOT.mkdir(parents=True, exist_ok=True)


def load_predictions(name: str) -> pd.DataFrame:
    pred_path = RESULTS_ROOT / f"{name}_predictions.tsv"
    if not pred_path.exists():
        raise FileNotFoundError(pred_path)
    df = pd.read_csv(pred_path, sep="\t")
    return df


def plot_roc(ax: plt.Axes, pred: pd.DataFrame, label: str, color: str) -> None:
    fpr, tpr, _ = roc_curve(pred["label"], pred["probability"])
    ax.plot(fpr, tpr, label=f"{label} (AUC={auc(fpr, tpr):.2f})", color=color, lw=2)


def figure_roc() -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    datasets = [
        ("msv000089062_csf", "CSF (MSV000089062)", sns.color_palette("tab10")[0]),
        ("pxd018301_csf", "CSF (PXD018301)", sns.color_palette("tab10")[2]),
        ("pxd032767_pxd018301_plasma", "Plasma (PXD032767 + PXD018301)", sns.color_palette("tab10")[1]),
    ]
    for name, label, color in datasets:
        pred = load_predictions(name)
        if pred["label"].nunique() < 2:
            continue
        plot_roc(ax, pred, label, color)

    ax.plot([0, 1], [0, 1], "--", color="grey", lw=1)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Liquid Site-Triage ROC Curves")
    ax.legend(loc="lower right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG_ROOT / "site_triage_roc.pdf")
    fig.savefig(FIG_ROOT / "site_triage_roc.png", dpi=300)
    plt.close(fig)


def figure_effects() -> None:
    robust_path = RESULTS_ROOT / "pxd032767_pxd018301_plasma_robust_logit.tsv"
    mixed_path = RESULTS_ROOT / "pxd032767_pxd018301_plasma_mixed_effects.tsv"

    robust = pd.read_csv(robust_path, sep="\t")
    robust = robust[robust["feature"] != "const"].copy()
    robust["type"] = "Robust Logit"
    robust["label"] = robust["feature"].replace({"dataset_indicator": "Dataset indicator"})

    mixed = pd.read_csv(mixed_path, sep="\t")
    mixed = mixed[mixed["feature"] != "Intercept"].copy()
    mixed["type"] = "Mixed Effects"
    mixed["label"] = mixed["feature"].replace({"dataset_indicator": "Dataset indicator"})

    df = pd.concat([
        robust.rename(columns={"coef": "estimate", "ci_lower": "lower", "ci_upper": "upper"})[
            ["label", "estimate", "lower", "upper", "type"]
        ],
        mixed.rename(columns={"coef": "estimate", "ci_lower": "lower", "ci_upper": "upper"})[
            ["label", "estimate", "lower", "upper", "type"]
        ],
    ])

    fig, ax = plt.subplots(figsize=(8, 5))
    sns.pointplot(
        data=df,
        y="label",
        x="estimate",
        hue="type",
        dodge=0.4,
        join=False,
        palette="tab10",
        ax=ax,
    )
    for _, row in df.iterrows():
        ax.plot([row["lower"], row["upper"]], [row["label"], row["label"]],
                color=sns.color_palette("tab10")[0 if row["type"] == "Robust Logit" else 1], lw=2)

    ax.axvline(0, color="grey", lw=1, linestyle="--")
    ax.set_xlabel("Coefficient (log-odds)")
    ax.set_ylabel("Feature")
    ax.set_title("Plasma Model Effect Sizes")
    ax.legend(title="Model", loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG_ROOT / "site_triage_effects.pdf")
    fig.savefig(FIG_ROOT / "site_triage_effects.png", dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    figure_roc()
    figure_effects()
    print("Figures written to", FIG_ROOT)


if __name__ == "__main__":
    main()
