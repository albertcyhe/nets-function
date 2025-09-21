# -*- coding: utf-8 -*-
"""Utilities for downloading named-metabolite tables from Metabolomics Workbench.

The main helper ``download_named_metabolites_via_api`` attempts to retrieve a
"named metabolite × sample" matrix via the official REST API. Whenever a
datatable is published for at least one analysis of a study, the function
returns a tidy ``pandas.DataFrame``. If no datatable exists, it optionally falls
back to reshaping the generic JSON ``/data`` endpoint.

Also included are convenience routines for harmonising reduced/oxidised
Glutathione labels (GSH / GSSG), which are useful when computing redox ratios.

References
---------
- REST documentation: https://www.metabolomicsworkbench.org/rest
- Example endpoints used:
    * ``/rest/study/study_id/<STXXXXXX>/analysis``
    * ``/rest/study/analysis_id/<ANXXXXXX>/datatable``
    * ``/rest/study/study_id/<STXXXXXX>/data``
"""

from __future__ import annotations

import io
import re
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd
import requests

__all__ = [
    "download_named_metabolites_via_api",
    "get_redox_alias_table",
    "normalize_name",
    "to_redox_canonical",
]


MW_BASE = "https://www.metabolomicsworkbench.org/rest"
HEADERS = {
    "User-Agent": "mw_parser/0.1 (+programmatic access for metabolomics preprocessing)",
    "Accept": "*/*",
}
STUDY_ID_RE = re.compile(r"^ST\d{6}$", re.IGNORECASE)
ANALYSIS_ID_RE = re.compile(r"^AN\d{6}$", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def download_named_metabolites_via_api(
    study_id: str,
    *,
    timeout: int = 60,
    session: Optional[requests.Session] = None,
    allow_fallback_to_data: bool = True,
) -> Optional[pd.DataFrame]:
    """Return a named-metabolite matrix for ``study_id`` if available via REST.

    Strategy
    --------
    1. List all ``analysis_id`` values published for the study.
    2. For each analysis, attempt ``/study/analysis_id/<AN>/datatable`` and
       parse it as a tab-delimited table. The first non-empty table is returned.
    3. If no datatable is available and ``allow_fallback_to_data`` is True,
       fetch ``/study/study_id/<ST>/data`` and reshape it into a matrix.
    4. If neither endpoint yields a usable matrix, return ``None``.
    """

    study_id = study_id.strip()
    if not STUDY_ID_RE.match(study_id):
        raise ValueError(f"Invalid study_id format: {study_id!r}")

    sess = session or requests.Session()

    analyses = _list_analyses_for_study(study_id, session=sess, timeout=timeout)
    for an_id in analyses:
        try:
            df = _fetch_datatable_for_analysis(an_id, session=sess, timeout=timeout)
        except Exception:
            df = None
        if _is_valid_matrix(df):
            df.attrs["study_id"] = study_id.upper()
            df.attrs["analysis_id"] = an_id.upper()
            df.attrs["source"] = f"{MW_BASE}/study/analysis_id/{an_id}/datatable"
            return df

    if allow_fallback_to_data:
        try:
            df = _fetch_study_data_and_pivot(study_id, session=sess, timeout=timeout)
        except Exception:
            df = None
        if _is_valid_matrix(df):
            df.attrs["study_id"] = study_id.upper()
            df.attrs["analysis_id"] = None
            df.attrs["source"] = f"{MW_BASE}/study/study_id/{study_id}/data"
            return df

    return None


# ---------------------------------------------------------------------------
# Internal helpers: REST handling
# ---------------------------------------------------------------------------


def _list_analyses_for_study(
    study_id: str,
    *,
    session: requests.Session,
    timeout: int,
) -> List[str]:
    url = f"{MW_BASE}/study/study_id/{study_id}/analysis"
    resp = session.get(url, headers=HEADERS, timeout=timeout)
    if resp.status_code != 200:
        return []

    try:
        payload = resp.json()
    except ValueError:
        return []

    raw_ids: List[str] = []
    if isinstance(payload, list):
        for row in payload:
            an_id = _extract_analysis_id(row)
            if an_id:
                raw_ids.append(an_id)
    elif isinstance(payload, dict):
        for value in payload.values():
            if isinstance(value, list):
                for row in value:
                    an_id = _extract_analysis_id(row)
                    if an_id:
                        raw_ids.append(an_id)

    # De-duplicate while preserving order
    ordered: List[str] = []
    seen = set()
    for an in raw_ids:
        if an and ANALYSIS_ID_RE.match(an) and an not in seen:
            seen.add(an)
            ordered.append(an)
    return ordered


def _extract_analysis_id(row: Dict[str, Any]) -> Optional[str]:
    for key in ("analysis_id", "ANALYSIS_ID", "Analysis_id", "analysisId"):
        if key in row and row[key]:
            return str(row[key]).strip()
    return None


def _fetch_datatable_for_analysis(
    analysis_id: str,
    *,
    session: requests.Session,
    timeout: int,
) -> Optional[pd.DataFrame]:
    if not ANALYSIS_ID_RE.match(analysis_id):
        return None
    url = f"{MW_BASE}/study/analysis_id/{analysis_id}/datatable"
    resp = session.get(url, headers=HEADERS, timeout=timeout)
    if resp.status_code != 200:
        return None

    text = resp.text or ""
    if not text.strip() or text.lower().startswith("error"):
        return None

    df = pd.read_csv(
        io.StringIO(text),
        sep="\t",
        header=0,
        dtype=str,
        engine="python",
        on_bad_lines="skip",
    )
    df = df.loc[:, [c for c in df.columns if str(c).strip() and not str(c).startswith("Unnamed")]]
    return df if not df.empty else None


def _fetch_study_data_and_pivot(
    study_id: str,
    *,
    session: requests.Session,
    timeout: int,
) -> Optional[pd.DataFrame]:
    url = f"{MW_BASE}/study/study_id/{study_id}/data"
    resp = session.get(url, headers=HEADERS, timeout=timeout)
    if resp.status_code != 200:
        return None

    try:
        payload = resp.json()
    except ValueError:
        return None

    records: Iterable[Dict[str, Any]]
    if isinstance(payload, list):
        records = payload
    elif isinstance(payload, dict):
        lists = [v for v in payload.values() if isinstance(v, list)]
        records = lists[0] if lists else []
    else:
        records = []

    rows: List[Tuple[str, str, Any]] = []
    meta_name_keys = ("metabolite_name", "METABOLITE_NAME", "name", "Name")
    data_keys = ("data", "DATA", "values", "VALUES", "results", "RESULTS")

    for rec in records:
        metabolite = None
        for key in meta_name_keys:
            if key in rec and rec[key]:
                metabolite = str(rec[key]).strip()
                break
        if not metabolite:
            for key in ("refmet_name", "REFMET_NAME"):
                if key in rec and rec[key]:
                    metabolite = str(rec[key]).strip()
                    break
        if not metabolite:
            continue

        container = None
        for key in data_keys:
            value = rec.get(key)
            if isinstance(value, dict):
                container = value
                break
        if container is None:
            for value in rec.values():
                if isinstance(value, dict) and _looks_like_sample_value_map(value):
                    container = value
                    break
        if not container:
            continue

        for sample_id, intensity in container.items():
            rows.append((metabolite, str(sample_id).strip(), intensity))

    if not rows:
        return None

    long_df = pd.DataFrame(rows, columns=["metabolite_name", "sample_id", "value"])
    wide = (
        long_df.pivot_table(index="metabolite_name", columns="sample_id", values="value", aggfunc="first")
        .sort_index()
        .reset_index()
    )
    return wide if not wide.empty else None


def _looks_like_sample_value_map(container: Dict[Any, Any]) -> bool:
    if not container:
        return False
    sample_keys = list(container.keys())[:10]
    sample_vals = [container[k] for k in sample_keys]
    keys_ok = all(isinstance(k, (str, int)) for k in sample_keys)
    vals_ok = all(
        (v is None)
        or isinstance(v, (int, float))
        or (isinstance(v, str) and bool(re.match(r"^\s*-?\d+(\.\d+)?(e[-+]?\d+)?\s*$", v)))
        for v in sample_vals
    )
    return keys_ok and vals_ok


def _is_valid_matrix(df: Optional[pd.DataFrame]) -> bool:
    return bool(df is not None and not df.empty and df.shape[1] >= 2)


# ---------------------------------------------------------------------------
# Redox alias mapping utilities (GSH / GSSG)
# ---------------------------------------------------------------------------


def normalize_name(name: str) -> str:
    """Lowercase + normalise punctuation/spelling."""
    if name is None:
        return ""
    name = name.strip().lower()
    name = name.replace("γ", "gamma")
    name = name.replace("oxidised", "oxidized")
    name = name.replace("disulphide", "disulfide")
    name = re.sub(r"[()\[\],;:_\-]+", " ", name)
    name = re.sub(r"\s+", " ", name)
    return name.strip()


REDOX_ALIAS_RECORDS: List[Dict[str, str]] = [
    # --- Reduced glutathione (GSH) ---
    {"alias": "GSH", "canonical": "GSH", "hmdb_id": "HMDB0000125", "kegg_id": "C00051", "chebi_id": "CHEBI:16856"},
    {"alias": "glutathione", "canonical": "GSH", "hmdb_id": "HMDB0000125", "kegg_id": "C00051", "chebi_id": "CHEBI:16856"},
    {"alias": "reduced glutathione", "canonical": "GSH", "hmdb_id": "HMDB0000125", "kegg_id": "C00051", "chebi_id": "CHEBI:177535"},
    {"alias": "glutathione reduced", "canonical": "GSH", "hmdb_id": "HMDB0000125", "kegg_id": "C00051", "chebi_id": "CHEBI:177535"},
    {"alias": "l-glutathione", "canonical": "GSH", "hmdb_id": "HMDB0000125", "kegg_id": "C00051", "chebi_id": "CHEBI:16856"},
    {"alias": "l-glutathione reduced", "canonical": "GSH", "hmdb_id": "HMDB0000125", "kegg_id": "C00051", "chebi_id": "CHEBI:177535"},
    {"alias": "l-glutathione (reduced)", "canonical": "GSH", "hmdb_id": "HMDB0000125", "kegg_id": "C00051", "chebi_id": "CHEBI:177535"},
    {"alias": "gamma-glutamyl-cysteinyl-glycine", "canonical": "GSH", "hmdb_id": "HMDB0000125", "kegg_id": "C00051", "chebi_id": "CHEBI:16856"},
    {"alias": "l-gamma-glutamyl-l-cysteinylglycine", "canonical": "GSH", "hmdb_id": "HMDB0000125", "kegg_id": "C00051", "chebi_id": "CHEBI:16856"},
    {"alias": "5-l-glutamyl-l-cysteinylglycine", "canonical": "GSH", "hmdb_id": "HMDB0000125", "kegg_id": "C00051", "chebi_id": "CHEBI:16856"},
    # --- Oxidized glutathione (GSSG) ---
    {"alias": "GSSG", "canonical": "GSSG", "hmdb_id": "HMDB0003337", "kegg_id": "C00127", "chebi_id": "CHEBI:17858"},
    {"alias": "oxidized glutathione", "canonical": "GSSG", "hmdb_id": "HMDB0003337", "kegg_id": "C00127", "chebi_id": "CHEBI:17858"},
    {"alias": "glutathione oxidized", "canonical": "GSSG", "hmdb_id": "HMDB0003337", "kegg_id": "C00127", "chebi_id": "CHEBI:17858"},
    {"alias": "glutathione disulfide", "canonical": "GSSG", "hmdb_id": "HMDB0003337", "kegg_id": "C00127", "chebi_id": "CHEBI:17858"},
    {"alias": "glutathione disulphide", "canonical": "GSSG", "hmdb_id": "HMDB0003337", "kegg_id": "C00127", "chebi_id": "CHEBI:17858"},
    {"alias": "oxiglutatione", "canonical": "GSSG", "hmdb_id": "HMDB0003337", "kegg_id": "C00127", "chebi_id": "CHEBI:17858"},
]

REDOX_ALIAS_MAP: Dict[str, str] = {normalize_name(rec["alias"]): rec["canonical"] for rec in REDOX_ALIAS_RECORDS}


def get_redox_alias_table() -> pd.DataFrame:
    """Return a DataFrame of redox aliases with canonical tags and IDs."""
    df = pd.DataFrame(REDOX_ALIAS_RECORDS)
    df["alias_norm"] = df["alias"].map(normalize_name)
    return df[["alias", "alias_norm", "canonical", "hmdb_id", "kegg_id", "chebi_id"]].copy()


def to_redox_canonical(name: str) -> Optional[str]:
    """Map *name* to canonical 'GSH'/'GSSG' when recognised, else ``None``."""
    if not name:
        return None
    return REDOX_ALIAS_MAP.get(normalize_name(name))
