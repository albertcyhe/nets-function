#!/usr/bin/env python3
"""Build liquid-site triage features from MaxQuant outputs or pre-quantified cohort tables.

Steps
-----
1. Parse a run manifest (`experiment`, `patient_id`, `bm_group`, ...)
2. Derive NE/PR3 footprint and THBS1/COL1A1/COL4A1 cleavage indices from `peptides.txt`
3. Summarise NET-F/NET-M/Serpin/THBS1 abundances from `proteinGroups.txt`
4. Merge into run-level feature table stored under `data/processed/proteomics/<dataset>/maxquant_run_features.tsv`
5. Optionally parse Supplementary Table S1-style Excel sheet to build patient-level features and train
   a Serpin + ΔFM logistic model (outputs saved to `results/tables`).

Usage (examples)
----------------
MaxQuant txt outputs + supplementary cohort table::

    micromamba run -n proteomics python scripts/py/site_triage_liquid.py \
      --mode maxquant \
      --dataset PXD032767 \
      --maxquant-dir data/interim/proteomics/PXD032767/txt \
      --run-manifest data/interim/proteomics/PXD032767/metadata/run_manifest.tsv \
      --clinical-xlsx data/interim/proteomics/PXD032767/metadata/vdac161_suppl_supplementary_table_s1.xlsx

CSF cohort XLSX (MSV000089062)::

    micromamba run -n proteomics python scripts/py/site_triage_liquid.py \
      --mode csf \
      --dataset MSV000089062 \
      --input-xlsx data/interim/proteomics/MSV000089062/metadata/vdac161_suppl_supplementary_table_s1.xlsx
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable, List, Optional, Sequence
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.calibration import calibration_curve
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    brier_score_loss,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
from statsmodels.stats.sandwich_covariance import cov_hc3

logger = logging.getLogger("site_triage_liquid")
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

MOTIF_P1 = set("AVILMF")
CTRL_P1 = set("DE")
SERPIN_CORE = ["SERPINA1", "SERPINA3", "SERPINI1"]
SERPIN_EXT = SERPIN_CORE + ["SLPI", "PI3", "ELAFIN", "A2M"]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["maxquant", "csf", "plasma"], default="maxquant", help="Processing mode")
    ap.add_argument("--dataset", required=True, help="Dataset identifier (e.g. PXD032767)")
    ap.add_argument("--maxquant-dir", type=Path, help="Path to MaxQuant txt/ directory")
    ap.add_argument("--run-manifest", type=Path, help="TSV mapping experiments to patients / groups")
    ap.add_argument(
        "--input-xlsx",
        type=Path,
        help="For --mode csf: Excel file containing clinical metadata and intensity matrices",
    )
    ap.add_argument(
        "--intensity-sheet",
        type=str,
        default="Raw Intensities",
        help="Sheet name holding the protein intensity matrix (CSF mode)",
    )
    ap.add_argument(
        "--scaled-sheet",
        type=str,
        default="Scaled Intensities",
        help="Optional sheet with scaled intensities (CSF mode)",
    )
    ap.add_argument(
        "--label-sheet",
        type=str,
        default="Clinical Metadata",
        help="Sheet name with clinical labels (CSF mode)",
    )
    ap.add_argument(
        "--label-tsv",
        type=Path,
        help="Optional TSV with sample annotations (columns: sample/patient_id, bm_group, ...)",
    )
    ap.add_argument(
        "--label-id-column",
        type=str,
        default="patient_id",
        help="Identifier column name inside --label-tsv (default: patient_id)",
    )
    ap.add_argument(
        "--primary-features",
        type=Path,
        help="Optional precomputed primary feature table (TSV) for plasma mode",
    )
    ap.add_argument(
        "--aux-dataset",
        type=str,
        help="Auxiliary dataset identifier for plasma mode",
    )
    ap.add_argument(
        "--aux-features",
        type=Path,
        help="Auxiliary feature table (TSV) for plasma mode",
    )
    ap.add_argument(
        "--aux-input-xlsx",
        type=Path,
        help="Optional auxiliary intensity matrix (if features need computing)",
    )
    ap.add_argument(
        "--aux-intensity-sheet",
        type=str,
        default="Raw Intensities",
        help="Sheet name for auxiliary intensity matrix",
    )
    ap.add_argument(
        "--aux-label-tsv",
        type=Path,
        help="Auxiliary label TSV",
    )
    ap.add_argument(
        "--clinical-xlsx",
        type=Path,
        default=None,
        help="Optional Supplementary Table S1 Excel for cohort-level intensities",
    )
    ap.add_argument(
        "--module-f",
        type=Path,
        default=Path("resources/modules/net_f_proteome.tsv"),
        help="NET-F proteomic module gene list (TSV with column 'Gene')",
    )
    ap.add_argument(
        "--module-m",
        type=Path,
        default=Path("resources/modules/net_m_proteome.tsv"),
        help="NET-M proteomic module gene list (TSV with column 'Gene')",
    )
    ap.add_argument(
        "--processed-root",
        type=Path,
        default=Path("data/processed/proteomics"),
        help="Directory to store processed feature tables",
    )
    ap.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results/tables"),
        help="Directory to store model summaries",
    )
    ap.add_argument("--folds", type=int, default=5, help="Number of CV folds for patient-level model")
    ap.add_argument("--random-state", type=int, default=42, help="Random seed for CV")
    return ap.parse_args()


def read_module_genes(path: Path) -> List[str]:
    df = pd.read_csv(path, sep="\t", header=None)
    genes = df.iloc[:, 0].astype(str).str.strip().tolist()
    return [g for g in genes if g]


def decision_curve(y_true: Sequence[int], probs: Sequence[float], thresholds: np.ndarray) -> pd.DataFrame:
    y_true = np.asarray(y_true)
    probs = np.asarray(probs)
    n = len(y_true)
    prevalence = y_true.mean()
    rows = []
    for thr in thresholds:
        if thr <= 0 or thr >= 1:
            continue
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


def merge_label_table(
    features: pd.DataFrame,
    label_path: Optional[Path],
    id_column: str = "patient_id",
) -> pd.DataFrame:
    if not label_path or not label_path.exists():
        return features
    label_df = pd.read_csv(label_path, sep="\t")
    if id_column not in label_df.columns:
        raise ValueError(f"Column '{id_column}' not found in label TSV {label_path}")
    label_df[id_column] = label_df[id_column].astype(str).str.strip()
    label_df = label_df.drop_duplicates(id_column)
    overlap = [c for c in label_df.columns if c != id_column and c in features.columns]
    if overlap:
        features = features.drop(columns=overlap)
    merged = features.join(label_df.set_index(id_column), how="left")
    logger.info("Merged manual label table %s", label_path)
    return merged


def load_manifest(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"experiment", "patient_id", "bm_group"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Run manifest missing columns: {sorted(missing)}")
    return df


def load_peptides(path: Path, experiments: Iterable[str]) -> pd.DataFrame:
    cols = [
        "Gene names",
        "Amino acid before",
        "First amino acid",
        "Last amino acid",
        "Amino acid after",
    ] + [f"Experiment {exp}" for exp in experiments]
    logger.info("Reading peptides from %s", path)
    df = pd.read_csv(path, sep="\t", usecols=lambda c: c in cols)
    for col in ["Gene names", "Amino acid before", "First amino acid", "Last amino acid", "Amino acid after"]:
        if col not in df.columns:
            df[col] = ""
    df[["Amino acid before", "First amino acid", "Last amino acid", "Amino acid after"]] = (
        df[["Amino acid before", "First amino acid", "Last amino acid", "Amino acid after"]]
        .fillna("")
        .astype(str)
        .apply(lambda s: s.str.strip().str.upper())
    )
    return df


def compute_semi_mask(df: pd.DataFrame) -> pd.Series:
    before = df["Amino acid before"]
    after = df["Amino acid after"]
    first = df["First amino acid"]
    last = df["Last amino acid"]

    is_tryptic_n = before.isin({"K", "R"}) & (first != "P")
    is_tryptic_n |= before.isin({"", "_", "-"})
    is_tryptic_c = after.isin({"K", "R"}) & (last != "P")
    is_tryptic_c |= after.isin({"", "_", "-"})

    n_enzy = is_tryptic_n.astype(int) + is_tryptic_c.astype(int)
    return n_enzy < 2


def compute_footprint(df: pd.DataFrame, experiments: Iterable[str]) -> pd.DataFrame:
    semi_mask = compute_semi_mask(df)
    motif_mask = semi_mask & df["First amino acid"].isin(MOTIF_P1) & (~df["Amino acid before"].isin({"K", "R"}))
    ctrl_mask = semi_mask & df["First amino acid"].isin(CTRL_P1) & (~df["Amino acid before"].isin({"K", "R"}))

    rows = []
    for exp in experiments:
        col = f"Experiment {exp}"
        if col not in df.columns:
            logger.warning("Experiment %s not present in peptides table", exp)
            continue
        detected = df[col].notna() & (df[col] != 0)
        semi = int((semi_mask & detected).sum())
        motif = int((motif_mask & detected).sum())
        ctrl = int((ctrl_mask & detected).sum())
        footprint = np.log2(motif + 1.0) - np.log2(semi + 1.0)
        ctrl_idx = np.log2(ctrl + 1.0) - np.log2(semi + 1.0)
        rows.append({
            "experiment": exp,
            "semi_count": semi,
            "motif_count": motif,
            "ctrl_count": ctrl,
            "footprint_index": footprint,
            "control_footprint_index": ctrl_idx,
        })
    return pd.DataFrame(rows)


def compute_thbs1(df: pd.DataFrame, experiments: Iterable[str]) -> pd.DataFrame:
    semi_mask = compute_semi_mask(df)

    gene_lists = (
        df["Gene names"].fillna("").astype(str).str.upper().str.split(";").apply(lambda items: [i.strip() for i in items if i])
    )

    def gene_mask(gene: str) -> pd.Series:
        gene_upper = gene.upper()
        return gene_lists.apply(lambda items: gene_upper in items)

    genes = ["THBS1", "COL1A1", "COL4A1"]
    records: List[dict] = []
    for exp in experiments:
        col = f"Experiment {exp}"
        if col not in df.columns:
            continue
        detected = df[col].notna() & (df[col] != 0)
        for gene in genes:
            mask = gene_mask(gene)
            coverage = int((mask & detected).sum())
            semi = int((mask & detected & semi_mask).sum())
            cleave = np.log2(semi + 1.0) - np.log2(coverage + 1.0)
            records.append({
                "experiment": exp,
                "gene": gene,
                "coverage": coverage,
                "semi": semi,
                "cleave_idx": cleave,
            })
    if not records:
        return pd.DataFrame(columns=["experiment"])

    tidy = pd.DataFrame(records)
    wide_idx = tidy.pivot(index="experiment", columns="gene", values="cleave_idx")
    wide_idx.columns = [f"{col.lower()}_cleave_idx" for col in wide_idx.columns]
    coverage = tidy.pivot(index="experiment", columns="gene", values="coverage")
    coverage.columns = [f"{col.lower()}_coverage" for col in coverage.columns]
    semi = tidy.pivot(index="experiment", columns="gene", values="semi")
    semi.columns = [f"{col.lower()}_semi" for col in semi.columns]

    out = pd.concat([wide_idx, coverage, semi], axis=1).reset_index()
    return out


def load_protein_groups(path: Path, experiments: Iterable[str]) -> pd.DataFrame:
    intensity_cols = {f"Intensity {exp}": exp for exp in experiments}
    lfq_cols = {f"LFQ intensity {exp}": exp for exp in experiments}
    usecols = ["Gene names", "Reverse", "Potential contaminant"] + list(intensity_cols.keys()) + list(lfq_cols.keys())
    logger.info("Reading proteinGroups from %s", path)
    df = pd.read_csv(path, sep="\t", usecols=lambda c: c in usecols)
    df = df[(df["Reverse"] != "+") & (df["Potential contaminant"] != "+")]
    df["Gene"] = df["Gene names"].fillna("").astype(str).str.split(";").str[0].str.strip()
    df = df[df["Gene"] != ""].copy()
    intensities = df[["Gene"] + list(intensity_cols.keys())]
    intensities = intensities.groupby("Gene").sum()
    intensities.columns = [intensity_cols[c] for c in intensities.columns]
    return intensities


def compute_module_features(intensity: pd.DataFrame, module_f: List[str], module_m: List[str]) -> pd.DataFrame:
    log_counts = np.log2(intensity + 1.0)
    present_f = [g for g in module_f if g in log_counts.index]
    present_m = [g for g in module_m if g in log_counts.index]
    netf = log_counts.loc[present_f].mean(axis=0) if present_f else pd.Series(np.nan, index=log_counts.columns)
    netm = log_counts.loc[present_m].mean(axis=0) if present_m else pd.Series(np.nan, index=log_counts.columns)

    def zscore(values: pd.Series) -> pd.Series:
        std = values.std(ddof=0)
        if std == 0 or np.isnan(std):
            return pd.Series(np.zeros_like(values), index=values.index)
        return (values - values.mean()) / std

    serpin_counts = intensity.reindex(SERPIN_EXT).fillna(0.0)
    clr = np.log(serpin_counts + 1.0) - np.log(serpin_counts + 1.0).mean(axis=0)
    core = clr.reindex(SERPIN_CORE).mean(axis=0)
    ext = clr.mean(axis=0)

    thbs1 = log_counts.loc["THBS1"] if "THBS1" in log_counts.index else pd.Series(np.nan, index=log_counts.columns)

    features = pd.DataFrame({
        "experiment": log_counts.columns,
        "netf_score_raw": netf.values,
        "netm_score_raw": netm.values,
        "netf_score_z": zscore(netf).values,
        "netm_score_z": zscore(netm).values,
        "proteo_deltafm": (zscore(netf) - zscore(netm)).values,
        "Serpin_score_core": core.values,
        "Serpin_score_ext": ext.values,
        "THBS1_log_intensity": thbs1.values,
    })
    features["n_netf_genes"] = len(present_f)
    features["n_netm_genes"] = len(present_m)
    return features


def build_run_features(args: argparse.Namespace) -> Path:
    manifest = load_manifest(args.run_manifest)
    experiments = manifest["experiment"].astype(str).tolist()

    peptides_path = args.maxquant_dir / "peptides.txt"
    protein_groups_path = args.maxquant_dir / "proteinGroups.txt"
    if not peptides_path.exists() or not protein_groups_path.exists():
        raise FileNotFoundError("Expected MaxQuant files peptides.txt / proteinGroups.txt not found")

    peptides = load_peptides(peptides_path, experiments)
    footprint = compute_footprint(peptides, experiments)
    cleavage = compute_thbs1(peptides, experiments)

    module_f = read_module_genes(args.module_f)
    module_m = read_module_genes(args.module_m)
    intensity = load_protein_groups(protein_groups_path, experiments)
    modules = compute_module_features(intensity, module_f, module_m)

    run_features = manifest.merge(footprint, on="experiment", how="left")
    run_features = run_features.merge(cleavage, on="experiment", how="left")
    run_features = run_features.merge(modules, on="experiment", how="left")

    out_dir = args.processed_root / args.dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "maxquant_run_features.tsv"
    run_features.to_csv(out_path, sep="\t", index=False)
    logger.info("Run-level features written to %s", out_path)
    return out_path


def load_clinical_table(path: Path) -> pd.DataFrame:
    xl = pd.ExcelFile(path)
    possible = ["Raw Intensities", "Scaled Intensities", "Clinical Metadata"]
    sheet = next((s for s in possible if s in xl.sheet_names), None)
    if sheet is None:
        raise ValueError(f"None of the expected sheets {possible} found in {path}")
    df = xl.parse(sheet)
    df = df.rename(columns={df.columns[0]: "Gene"})
    df["Gene"] = df["Gene"].astype(str).str.strip()
    df = df[df["Gene"] != ""]
    df = df.set_index("Gene")
    return df


def build_patient_features(args: argparse.Namespace) -> Optional[pd.DataFrame]:
    if args.clinical_xlsx is None or not args.clinical_xlsx.exists():
        logger.info("Clinical XLSX not provided; skipping patient-level processing")
        return None

    logger.info("Loading cohort matrix from %s", args.clinical_xlsx)
    mat = load_clinical_table(args.clinical_xlsx)
    mat = mat.clip(lower=1.0)

    module_f = read_module_genes(args.module_f)
    module_m = read_module_genes(args.module_m)

    log_mat = np.log(mat + 1.0)
    netf = log_mat.loc[[g for g in module_f if g in log_mat.index]].mean(axis=0)
    netm = log_mat.loc[[g for g in module_m if g in log_mat.index]].mean(axis=0)

    def zscore(series: pd.Series) -> pd.Series:
        std = series.std(ddof=0)
        if std == 0 or np.isnan(std):
            return pd.Series(np.zeros_like(series), index=series.index)
        return (series - series.mean()) / std

    serpin_genes = [g for g in SERPIN_CORE if g in mat.index]
    if serpin_genes:
        clr = np.log(mat.loc[serpin_genes] + 1.0) - np.log(mat.loc[serpin_genes] + 1.0).mean(axis=0)
        serpin_core = clr.mean(axis=0)
    else:
        serpin_core = pd.Series(np.nan, index=mat.columns)

    deltafm = zscore(netf) - zscore(netm)

    patient_features = pd.DataFrame({
        "patient_id": mat.columns,
        "Serpin_score_core": serpin_core.values,
        "Proteo_DeltaFM": deltafm.values,
    }).set_index("patient_id")
    return patient_features


def attach_labels(patient_features: pd.DataFrame, clinical_sheet: Path) -> pd.DataFrame:
    xl = pd.ExcelFile(clinical_sheet)
    if "Clinical Metadata" not in xl.sheet_names:
        logger.warning("Clinical Metadata sheet not found; labels unavailable")
        return patient_features
    clin = xl.parse("Clinical Metadata")
    clin.columns = clin.columns.str.strip().str.lower().str.replace(" ", "_")
    clin = clin.rename(columns={"patient_id": "patient_id"})
    clin["patient_id"] = clin["patient_id"].astype(str).str.strip()
    if "diagnosis" in clin.columns:
        clin["diagnosis"] = clin["diagnosis"].astype(str).str.strip()
    if "bm_group" not in clin.columns:
        def assign(row):
            diag = str(row.get("diagnosis", "")).lower()
            if "brainmet" in diag or diag in {"brainm", "brain"}:
                return "brain_met"
            if diag in {"healthy", "control", "nph"}:
                return "non_brain_control"
            return "non_brain_other"
        clin["bm_group"] = clin.apply(assign, axis=1)
    clin = clin[["patient_id", "bm_group"]]
    clin = clin.dropna(subset=["patient_id"]).drop_duplicates()
    clin = clin.set_index("patient_id")
    merged = patient_features.join(clin, how="left")
    return merged


def build_csf_features(
    args: argparse.Namespace,
    module_f: List[str],
    module_m: List[str],
) -> pd.DataFrame:
    if args.input_xlsx is None or not args.input_xlsx.exists():
        raise FileNotFoundError("--input-xlsx is required for CSF mode")

    xl = pd.ExcelFile(args.input_xlsx)
    multi_header = False
    try:
        df = xl.parse(args.intensity_sheet, header=[0, 1, 2])
        if isinstance(df.columns, pd.MultiIndex):
            cols = df.columns
            area_cols_check = [col for col in cols if "area" in str(col[-1]).lower()]
            if area_cols_check:
                multi_header = True
            else:
                flattened = []
                for col in cols:
                    levels = [str(level).strip() for level in col if str(level).strip() and not str(level).lower().startswith("unnamed")]
                    flattened.append(levels[-1] if levels else str(col[-1]))
                df.columns = flattened
        else:
            df = pd.DataFrame()
    except ValueError:
        df = pd.DataFrame()

    if not multi_header or df.empty:
        df = xl.parse(args.intensity_sheet, header=0)
        if not df.empty:
            first_val = str(df.iloc[0, 0]).strip().lower() if pd.notna(df.iloc[0, 0]) else ""
            if first_val in {"", "source:"}:
                df = df.iloc[1:, :]
            if not df.empty:
                first_val = str(df.iloc[0, 0]).strip().lower()
                if first_val in {"accession", "gene"}:
                    df = df.iloc[1:, :]
        df = df.reset_index(drop=True)

    if isinstance(df.columns, pd.MultiIndex):
        cols = df.columns

        def find_col(name: str):
            name = name.lower()
            for col in cols:
                if str(col[2]).strip().lower() == name:
                    return col
            return None

        acc_col = find_col("accession") or cols[0]
        desc_col = find_col("description")
        area_cols = [
            col for col in cols if "area" in str(col[2]).lower()
        ]
        if not area_cols:
            raise ValueError("No Area columns detected in multi-index header")

        accession = df[acc_col].astype(str)
        description = df[desc_col].astype(str) if desc_col else pd.Series(index=df.index, dtype=str)
        gene_symbol = description.str.extract(r"GN=([A-Za-z0-9\-]+)")[0]
        genes = gene_symbol.fillna(accession).str.strip()
        genes = genes[genes != ""]

        intensity = df[area_cols].copy()
        sample_names = []
        used = {}
        for col in area_cols:
            candidate = str(col[1]).strip()
            if not candidate or candidate.lower().startswith("unnamed"):
                candidate = str(col[2]).split(":")[0].strip()
            candidate = candidate or "sample"
            base = candidate
            idx = 1
            while candidate in used:
                idx += 1
                candidate = f"{base}_{idx}"
            used[candidate] = True
            sample_names.append(candidate)
        intensity.columns = sample_names
        intensity = intensity.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        intensity.index = genes
        intensity = intensity[intensity.index != ""]
        intensity = intensity[~intensity.index.duplicated(keep="first")]
    else:
        columns = list(df.columns)
        if columns:
            df = df.rename(columns={columns[0]: "Accession"})
        if len(columns) > 1:
            df = df.rename(columns={columns[1]: "Description"})

        if "Description" in df.columns:
            gene_symbol = df["Description"].astype(str).str.extract(r"GN=([A-Za-z0-9\-]+)")[0]
        else:
            gene_symbol = pd.Series(dtype=str)

        df["Gene"] = gene_symbol.fillna(df.get("Accession", "")).astype(str).str.strip()
        df = df[df["Gene"] != ""]
        df = df[~df["Gene"].str.lower().isin({"accession", "gene"})]

        numeric_df = df.drop(columns=["Description"], errors="ignore")
        intensity = numeric_df.set_index("Gene").apply(pd.to_numeric, errors="coerce").fillna(0.0)
        logger.info("Initial columns detected: %s", intensity.columns[:5].tolist())
        valid_columns = [
            c
            for c in intensity.columns
            if isinstance(c, str)
            and c.strip() != ""
            and not c.lower().startswith("unnamed")
            and not c.lower().startswith("source")
            and c.lower() not in {"accession", "gene"}
        ]
        intensity = intensity[valid_columns]
        intensity = intensity[~intensity.index.duplicated(keep="first")]

    log_intensity = np.log2(intensity + 1.0)
    logger.info("log_intensity shape %s", log_intensity.shape)
    if log_intensity.shape[1] > 0:
        logger.info("Sample name preview: %s", list(log_intensity.columns[:5]))

    module_f_present = [g for g in module_f if g in log_intensity.index]
    module_m_present = [g for g in module_m if g in log_intensity.index]
    logger.info("Module genes present (F=%d, M=%d)", len(module_f_present), len(module_m_present))

    def zscore(series: pd.Series) -> pd.Series:
        std = series.std(ddof=0)
        if std == 0 or np.isnan(std):
            return pd.Series(np.zeros_like(series), index=series.index)
        return (series - series.mean()) / std

    netf_raw = log_intensity.loc[module_f_present].mean(axis=0) if module_f_present else pd.Series(np.nan, index=log_intensity.columns)
    netm_raw = log_intensity.loc[module_m_present].mean(axis=0) if module_m_present else pd.Series(np.nan, index=log_intensity.columns)
    netf_z = zscore(netf_raw.dropna()).reindex(log_intensity.columns)
    netm_z = zscore(netm_raw.dropna()).reindex(log_intensity.columns)

    serpin_genes = [g for g in SERPIN_EXT if g in intensity.index]
    logger.info("Serpin genes present: %d", len(serpin_genes))
    if serpin_genes:
        serpin_counts = intensity.loc[serpin_genes]
        log_serpin = np.log(serpin_counts + 1.0)
        clr = log_serpin.subtract(log_serpin.mean(axis=0), axis=1)
        serpin_core = clr.loc[[g for g in SERPIN_CORE if g in clr.index]].mean(axis=0)
        serpin_ext = clr.mean(axis=0)
    else:
        serpin_core = pd.Series(np.nan, index=log_intensity.columns)
        serpin_ext = pd.Series(np.nan, index=log_intensity.columns)

    def safe_mean(genes: Sequence[str]) -> pd.Series:
        present = [g for g in genes if g in log_intensity.index]
        if not present:
            return pd.Series(np.nan, index=log_intensity.columns)
        return log_intensity.loc[present].mean(axis=0)

    footprint_pos = safe_mean(["ELANE", "PRTN3", "CTSG", "MPO"])
    footprint_ctrl = safe_mean(["LAMB1", "LAMB2", "LAMA2", "LAMC1", "ELN"])
    footprint_index = footprint_pos - footprint_ctrl

    thbs1 = log_intensity.loc["THBS1"] if "THBS1" in log_intensity.index else pd.Series(np.nan, index=log_intensity.columns)
    collagen_mean = safe_mean(["COL1A1", "COL4A1", "COL4A2", "COL4A5"])
    thbs1_cleave = thbs1 - collagen_mean

    features = pd.DataFrame({
        "Serpin_score_core": serpin_core,
        "Serpin_score_ext": serpin_ext,
        "NetF_score_z": netf_z,
        "NetM_score_z": netm_z,
        "Proteo_DeltaFM": netf_z - netm_z,
        "THBS1_log": thbs1,
        "footprint_index": footprint_index,
        "THBS1_cleave_idx": thbs1_cleave,
    })

    for col in ["THBS1_log", "footprint_index", "THBS1_cleave_idx"]:
        if col in features.columns:
            if features[col].isna().all():
                features[col] = 0.0
            else:
                features[col] = features[col].fillna(0.0)

    features = features.fillna(0.0)
    features["sample_id"] = features.index.astype(str)
    if "dataset" not in features.columns:
        features["dataset"] = args.dataset
    if "matrix" not in features.columns:
        features["matrix"] = "csf"
    logger.info("Raw feature matrix shape before labels %s", features.shape)

    if args.label_sheet in xl.sheet_names:
        clinical = xl.parse(args.label_sheet)
        clinical.columns = clinical.columns.str.strip().str.lower().str.replace(" ", "_")
        id_col = "patient_id" if "patient_id" in clinical.columns else clinical.columns[0]
        clinical["patient_id"] = clinical[id_col].astype(str).str.strip()

        def assign(row):
            diagnosis = str(row.get("diagnosis", "")).strip().lower()
            mapping = {
                "brainmet": "brain_met",
                "brain metastasis": "brain_met",
                "bm": "brain_met",
                "gbm": "non_brain_gbm",
                "pcnsl": "non_brain_cnsl",
                "cns lymphoma": "non_brain_cnsl",
                "nph": "non_brain_control",
            }
            return mapping.get(diagnosis, "non_brain_other")

        clinical["bm_group"] = clinical.apply(assign, axis=1)
        clinical = clinical.set_index("patient_id")
        features = features.join(clinical, how="left")
    else:
        logger.warning("Label sheet %s not found in %s", args.label_sheet, args.input_xlsx)

    logger.info("Feature matrix shape before manual merge %s", features.shape)
    features = merge_label_table(features, args.label_tsv, args.label_id_column)

    out_dir = args.processed_root / args.dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "csf_patient_features.tsv"
    logger.info("CSF feature matrix shape %s", features.shape)
    features.to_csv(out_path, sep="\t")
    logger.info("CSF patient-level features written to %s", out_path)
    return features


def train_patient_model(
    features: pd.DataFrame,
    feature_cols: Sequence[str],
    folds: int,
    random_state: int,
    group_col: Optional[str] = None,
) -> Optional[dict]:
    """Fit a logistic model with cross-validated probabilities.

    Returns a dictionary with summary statistics, cross-validated probabilities,
    fitted coefficients, and scaler parameters, or ``None`` if training is not
    feasible (e.g., labels missing).
    """

    if features is None or features.empty:
        return None
    if "bm_group" not in features.columns:
        logger.warning("bm_group labels missing; skipping model training")
        return None

    available_cols = [col for col in feature_cols if col in features.columns]
    if not available_cols:
        logger.warning("No requested feature columns available; skipping model training")
        return None

    labelled = features.dropna(subset=["bm_group"] + available_cols).copy()
    if labelled.empty:
        logger.warning("No labelled samples after filtering; skipping model training")
        return None

    labelled["label"] = (labelled["bm_group"].str.lower() == "brain_met").astype(int)
    if labelled["label"].nunique() < 2:
        logger.warning("Need at least one positive and negative sample; skipping model training")
        return None

    X = labelled[available_cols].astype(float).values
    y = labelled["label"].values
    indices = labelled.index.to_numpy()
    groups = None
    if group_col and group_col in features.columns:
        group_values = features.loc[indices, group_col].astype(str).fillna("unknown")
        if group_values.nunique() >= 2:
            groups = group_values

    n_splits = min(folds, labelled["label"].value_counts().min())
    if n_splits < 2:
        logger.warning("Insufficient samples per class for cross-validation; skipping model training")
        return None

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    cv_probs = np.zeros_like(y, dtype=float)
    cv_preds = np.zeros_like(y, dtype=int)
    aucs: List[float] = []
    accs: List[float] = []
    pr_aucs: List[float] = []
    briers: List[float] = []

    for train_idx, test_idx in cv.split(X, y):
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train_idx])
        X_test = scaler.transform(X[test_idx])
        model = LogisticRegression(max_iter=1000)
        model.fit(X_train, y[train_idx])
        probs = model.predict_proba(X_test)[:, 1]
        preds = model.predict(X_test)
        cv_probs[test_idx] = probs
        cv_preds[test_idx] = preds
        aucs.append(roc_auc_score(y[test_idx], probs))
        accs.append(accuracy_score(y[test_idx], preds))
        pr_aucs.append(average_precision_score(y[test_idx], probs))
        briers.append(brier_score_loss(y[test_idx], probs))

    final_scaler = StandardScaler().fit(X)
    X_scaled = final_scaler.transform(X)
    final_model = LogisticRegression(max_iter=1000)
    final_model.fit(X_scaled, y)
    scaled_df = pd.DataFrame(X_scaled, columns=available_cols, index=indices)

    summary = pd.DataFrame({
        "n_samples": [len(y)],
        "n_positive": [int(y.sum())],
        "n_negative": [int((y == 0).sum())],
        "cv_auc_mean": [float(np.mean(aucs))],
        "cv_auc_sd": [float(np.std(aucs))],
        "cv_accuracy_mean": [float(np.mean(accs))],
        "cv_accuracy_sd": [float(np.std(accs))],
        "cv_pr_auc_mean": [float(np.mean(pr_aucs))],
        "cv_pr_auc_sd": [float(np.std(pr_aucs))],
        "cv_brier_mean": [float(np.mean(briers))],
        "cv_brier_sd": [float(np.std(briers))],
    })

    mixed_info = None
    if groups is not None:
        try:
            data = scaled_df.copy()
            data["label"] = y
            data["group"] = groups.values
            formula = "label ~ " + " + ".join(available_cols)
            vc_formula = {"group": "0 + C(group)"}
            mixed_model = sm.BinomialBayesMixedGLM.from_formula(formula, vc_formula, data)
            mixed_fit = mixed_model.fit_map()
            fe_names = mixed_fit.model.exog_names
            fe_params = pd.Series(mixed_fit.fe_mean, index=fe_names)
            fe_se = pd.Series(mixed_fit.fe_sd, index=fe_names)
            mixed_info = {
                "params": fe_params,
                "bse": fe_se,
            }
        except Exception as exc:
            logger.warning("Mixed-effects model failed: %s", exc)

    return {
        "summary": summary,
        "y_true": pd.Series(y, index=indices, name="label"),
        "cv_probs": pd.Series(cv_probs, index=indices, name="prob"),
        "cv_preds": pd.Series(cv_preds, index=indices, name="pred"),
        "feature_names": available_cols,
        "coefficients": pd.Series(final_model.coef_[0], index=available_cols, name="coefficient"),
        "intercept": float(final_model.intercept_[0]),
        "scaler_mean": pd.Series(final_scaler.mean_, index=available_cols, name="mean"),
        "scaler_scale": pd.Series(final_scaler.scale_, index=available_cols, name="scale"),
        "mixed": mixed_info,
        "group_col": group_col,
    }


def export_model_outputs(
    dataset: str,
    suffix: str,
    result: dict,
    features: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    args.results_dir.mkdir(parents=True, exist_ok=True)
    base = f"{dataset.lower()}_{suffix}"

    summary = result["summary"].copy()
    y_true = result["y_true"].astype(int)
    probs = result["cv_probs"].astype(float)
    preds = result["cv_preds"].astype(int)

    auc_macro = roc_auc_score(y_true, probs)
    ap_macro = average_precision_score(y_true, probs)
    brier_macro = brier_score_loss(y_true, probs)
    summary["auc_overall"] = auc_macro
    summary["ap_overall"] = ap_macro
    summary["brier_overall"] = brier_macro

    summary_path = args.results_dir / f"{base}_metrics.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    logger.info("Model summary written to %s", summary_path)

    coef_df = result["coefficients"].reset_index()
    coef_df.columns = ["feature", "coefficient"]
    intercept_row = pd.DataFrame({"feature": ["intercept"], "coefficient": [result["intercept"]]})
    coef_df = pd.concat([coef_df, intercept_row], ignore_index=True)
    coef_path = args.results_dir / f"{base}_coefficients.tsv"
    coef_df.to_csv(coef_path, sep="\t", index=False)

    scaler = pd.concat([result["scaler_mean"], result["scaler_scale"]], axis=1)
    scaler.columns = ["mean", "scale"]
    scaler_path = args.results_dir / f"{base}_scaler.tsv"
    scaler.to_csv(scaler_path, sep="\t")

    preds_df = pd.DataFrame({
        "patient_id": y_true.index,
        "label": y_true.values,
        "probability": probs.values,
        "prediction": preds.values,
    })
    if "bm_group" in features.columns:
        preds_df["bm_group"] = features.reindex(y_true.index)["bm_group"].values
    preds_path = args.results_dir / f"{base}_predictions.tsv"
    preds_df.to_csv(preds_path, sep="\t", index=False)

    fpr, tpr, roc_thresh = roc_curve(y_true, probs)
    roc_df = pd.DataFrame({"threshold": roc_thresh, "fpr": fpr, "tpr": tpr})
    roc_df.to_csv(args.results_dir / f"{base}_roc.tsv", sep="\t", index=False)

    precision, recall, pr_thresh = precision_recall_curve(y_true, probs)
    pr_df = pd.DataFrame({"recall": recall, "precision": precision})
    pr_df.to_csv(args.results_dir / f"{base}_pr.tsv", sep="\t", index=False)

    try:
        prob_true, prob_pred = calibration_curve(y_true, probs, n_bins=min(10, max(3, int(len(y_true) / 5))), strategy="quantile")
        calib_df = pd.DataFrame({"prob_true": prob_true, "prob_pred": prob_pred})
        calib_df.to_csv(args.results_dir / f"{base}_calibration.tsv", sep="\t", index=False)
    except Exception as exc:
        logger.warning("Calibration curve failed: %s", exc)

    thresholds = np.linspace(0.01, 0.99, 99)
    dca_df = decision_curve(y_true.values, probs.values, thresholds)
    dca_df.to_csv(args.results_dir / f"{base}_dca.tsv", sep="\t", index=False)

    mixed = result.get("mixed")
    if mixed:
        params = mixed["params"]
        se = mixed["bse"]
        z = params / se
        p = 2 * stats.norm.sf(np.abs(z))
        ci_lower = params - 1.96 * se
        ci_upper = params + 1.96 * se
        coef_arr = np.asarray(params)
        se_arr = np.asarray(se)
        z_arr = np.asarray(z)
        p_arr = np.asarray(p)
        ci_lower_arr = np.asarray(ci_lower)
        ci_upper_arr = np.asarray(ci_upper)
        mixed_df = pd.DataFrame({
            "feature": params.index,
            "coef": coef_arr,
            "se": se_arr,
            "z": z_arr,
            "p_value": p_arr,
            "ci_lower": ci_lower_arr,
            "ci_upper": ci_upper_arr,
        })
        mixed_path = args.results_dir / f"{base}_mixed_effects.tsv"
        mixed_df.to_csv(mixed_path, sep="\t", index=False)


def export_association_tables(
    dataset: str,
    suffix: str,
    features: pd.DataFrame,
    feature_cols: Sequence[str],
    results_dir: Path,
) -> None:
    available_cols = [col for col in feature_cols if col in features.columns]
    if not available_cols:
        return

    df = features.dropna(subset=["bm_group"] + available_cols).copy()
    if df.empty:
        return

    df["label"] = (df["bm_group"].astype(str).str.lower() == "brain_met").astype(int)

    base = f"{dataset.lower()}_{suffix}"
    results_dir.mkdir(parents=True, exist_ok=True)

    spearman_rows = []
    for col in available_cols:
        if df[col].nunique() <= 1:
            continue
        rho, pval = stats.spearmanr(df[col], df["label"])
        spearman_rows.append({
            "feature": col,
            "spearman_rho": rho,
            "spearman_p": pval,
            "n": len(df),
        })
    if spearman_rows:
        spearman_df = pd.DataFrame(spearman_rows)
        spearman_df.to_csv(results_dir / f"{base}_spearman.tsv", sep="\t", index=False)

    X = sm.add_constant(df[available_cols], has_constant="add")
    try:
        logit_model = sm.Logit(df["label"], X)
        res = logit_model.fit(disp=0)
        mu = res.predict()
        y_array = df["label"].to_numpy()
        resid = y_array - mu
        w = mu * (1 - mu)
        XtWX_inv = np.asarray(res.normalized_cov_params)
        X_np = X.to_numpy()
        h = np.einsum('ij,jk,ik->i', X_np, XtWX_inv, X_np) * w
        wt = (resid / (1 - h)) ** 2
        wt_np = np.asarray(wt)
        meat = X_np.T @ (wt_np[:, None] * X_np)
        robust_cov = XtWX_inv @ meat @ XtWX_inv
        params = res.params
        se = np.sqrt(np.diag(robust_cov))
        z = params / se
        p_vals = 2 * stats.norm.sf(np.abs(z))
        ci_lower = params - 1.96 * se
        ci_upper = params + 1.96 * se
        robust_df = pd.DataFrame({
            "feature": params.index,
            "coef": params.values,
            "robust_se": se,
            "z": z,
            "p_value": p_vals,
            "ci_lower": ci_lower,
            "ci_upper": ci_upper,
        })
        robust_df.to_csv(results_dir / f"{base}_robust_logit.tsv", sep="\t", index=False)
    except Exception as exc:
        logger.warning("Robust logistic regression failed: %s", exc)


def load_feature_table(path: Path, dataset: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Feature table not found: {path}")
    df = pd.read_csv(path, sep="\t", index_col=0)
    df["sample_id"] = df.index.astype(str)
    df["dataset"] = dataset
    return df


def prepare_plasma_features(df: pd.DataFrame, dataset: str) -> pd.DataFrame:
    rename_map = {
        "THBS1_log_intensity": "THBS1_log",
        "THBS1_log": "THBS1_log",
    }
    df = df.rename(columns=rename_map)

    required = [
        "Serpin_score_core",
        "Proteo_DeltaFM",
        "THBS1_log",
        "footprint_index",
        "THBS1_cleave_idx",
    ]
    for col in required:
        if col not in df.columns:
            df[col] = 0.0

    if "bm_group" not in df.columns:
        df["bm_group"] = np.nan

    if "matrix" not in df.columns:
        df["matrix"] = "plasma_ev"

    df["dataset"] = dataset
    if "sample_id" not in df.columns:
        df["sample_id"] = df.index.astype(str)
    return df


def main() -> None:
    args = parse_args()

    module_f = read_module_genes(args.module_f)
    module_m = read_module_genes(args.module_m)

    if args.mode == "maxquant":
        if args.maxquant_dir is None or args.run_manifest is None:
            raise ValueError("--maxquant-dir and --run-manifest are required for maxquant mode")

        run_features_path = build_run_features(args)

        patient_features = build_patient_features(args)
        if patient_features is not None:
            patient_features = attach_labels(patient_features, args.clinical_xlsx)
            patient_features = merge_label_table(patient_features, args.label_tsv, args.label_id_column)
            out_dir = args.processed_root / args.dataset
            patient_path = out_dir / "maxquant_patient_features.tsv"
            patient_features.to_csv(patient_path, sep="\t")
            logger.info("Patient-level features written to %s", patient_path)

            feature_cols = [
                "Serpin_score_core",
                "Proteo_DeltaFM",
                "THBS1_log_intensity",
                "footprint_index",
                "THBS1_cleave_idx",
            ]
            result = train_patient_model(patient_features, feature_cols, args.folds, args.random_state, group_col=None)
            if result is not None:
                export_model_outputs(args.dataset, "liquid", result, patient_features, args)
                export_association_tables(args.dataset, "liquid", patient_features, feature_cols, args.results_dir)
            else:
                logger.info("Skipping model summary (insufficient label diversity)")

        logger.info("Done. Run features: %s", run_features_path)

    elif args.mode == "csf":
        features = build_csf_features(args, module_f, module_m)
        feature_cols = [
            "Serpin_score_core",
            "THBS1_log",
            "footprint_index",
            "THBS1_cleave_idx",
            "Proteo_DeltaFM",
        ]
        result = train_patient_model(features, feature_cols, args.folds, args.random_state, group_col=None)
        if result is not None:
            export_model_outputs(args.dataset, "csf", result, features, args)
            export_association_tables(args.dataset, "csf", features, feature_cols, args.results_dir)
        else:
            logger.info("Skipping CSF model summary (insufficient label diversity)")

    elif args.mode == "plasma":
        if not args.aux_dataset:
            raise ValueError("--aux-dataset is required for plasma mode")

        primary_path = args.primary_features or (args.processed_root / args.dataset / "maxquant_patient_features.tsv")
        if not primary_path.exists():
            logger.info("Primary feature table %s not found; computing via maxquant mode", primary_path)
            if args.maxquant_dir is None or args.run_manifest is None:
                raise FileNotFoundError("Primary features missing and maxquant inputs not provided")
            run_features_path = build_run_features(args)
            patient_features = build_patient_features(args)
            if patient_features is None:
                raise RuntimeError("Unable to derive primary patient features")
            patient_features = attach_labels(patient_features, args.clinical_xlsx)
            patient_features = merge_label_table(patient_features, args.label_tsv, args.label_id_column)
            out_dir = args.processed_root / args.dataset
            out_dir.mkdir(parents=True, exist_ok=True)
            patient_path = out_dir / "maxquant_patient_features.tsv"
            patient_features.to_csv(patient_path, sep="\t")
            primary_path = patient_path
            logger.info("Primary patient-level features written to %s", primary_path)

        primary_df = load_feature_table(primary_path, args.dataset)
        primary_df = merge_label_table(primary_df, args.label_tsv, args.label_id_column)
        primary_df = prepare_plasma_features(primary_df, args.dataset)

        aux_path = args.aux_features or (args.processed_root / args.aux_dataset / "csf_patient_features.tsv")
        if not aux_path.exists():
            if args.aux_input_xlsx is None:
                raise FileNotFoundError(
                    f"Auxiliary features not found ({aux_path}). Provide --aux-input-xlsx or precompute the table."
                )
            aux_args = SimpleNamespace(
                input_xlsx=args.aux_input_xlsx,
                intensity_sheet=args.aux_intensity_sheet,
                label_sheet=args.label_sheet,
                label_tsv=args.aux_label_tsv,
                label_id_column=args.label_id_column,
                processed_root=args.processed_root,
                dataset=args.aux_dataset,
            )
            aux_features = build_csf_features(aux_args, module_f, module_m)
            aux_path.parent.mkdir(parents=True, exist_ok=True)
            aux_features.to_csv(aux_path, sep="\t")
            logger.info("Auxiliary features written to %s", aux_path)

        aux_df = load_feature_table(aux_path, args.aux_dataset)
        aux_df = merge_label_table(aux_df, args.aux_label_tsv, args.label_id_column)
        aux_df = prepare_plasma_features(aux_df, args.aux_dataset)

        combined = pd.concat([primary_df, aux_df], axis=0, ignore_index=False)
        combined = combined.dropna(subset=["bm_group"])
        combined["dataset_indicator"] = (combined["dataset"].astype(str) == args.aux_dataset).astype(int)

        feature_cols = [
            "Serpin_score_core",
            "Proteo_DeltaFM",
            "THBS1_log",
            "footprint_index",
            "THBS1_cleave_idx",
            "dataset_indicator",
        ]

        result = train_patient_model(combined, feature_cols, args.folds, args.random_state, group_col="dataset")
        if result is not None:
            combined_name = f"{args.dataset}_{args.aux_dataset}"
            export_model_outputs(combined_name, "plasma", result, combined, args)
            export_association_tables(combined_name, "plasma", combined, feature_cols, args.results_dir)
        else:
            logger.info("Skipping plasma model summary (insufficient label diversity)")

    else:
        raise ValueError(f"Unsupported mode: {args.mode}")


if __name__ == "__main__":
    main()
