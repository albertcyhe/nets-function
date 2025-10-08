Data Sources and One-Click Retrieval
====================================

This document catalogs all external datasets and reference resources used across the project and provides a single entry point to re‑materialize them on a new server.

Quick start (one-click):
- Run `bash scripts/sh/download_all_data.sh` to pull the core datasets and references.
- See `README.md` for how to use these files in the analysis steps.

Proteomics (PRIDE / MassIVE)
- PRIDE: `PXD011796` — NET reference proteome (label-free). Used for module definition and validation.
- PRIDE: `PXD005719` — Brain metastasis membrane proteome; BBB crossing cargo focus.
- PRIDE: `PXD046330` — HER2+ model exposed to brain conditioned media (secretome).
- PRIDE: `PXD051579` — HER2+ brain metastasis BBB dysfunction (paused in cross-omics until matching transcriptome exists).
- PRIDE: `PXD032767` — Plasma proteomics used in the site-triage prototype (non-brain RAW not publicly exposed; use Supplementary Table S1 for labels and MaxQuant export for features).
- PRIDE: `PXD018301` — Plasma EV proteomics, used as a non-brain reference in triage fusion.
- MassIVE: `MSV000089062` — CSF cohort for triage analyses (summary tables are used in results; RAW/ID results retrieval is manual). 

Transcriptomics (GEO / Xena)
- GEO bulk paired cohorts: `GSE184869` (RNA-seq; processed log2 TMM CPM used), `GSE125989` (Affymetrix GPL571).
- GEO cohorts for metastasis stratification: `GSE43837` (brain), `GSE14017` and `GSE14018` (non-brain).
- GEO bigWig tracks: `GSE96860` (SKBR3 references for cross-omics via bigWig summarization).
- GEO single-cell: `GSE186344` (supplementary downloads used; RAW tar was unreliable).
- UCSC Xena: `MET500` expression and clinical annotations (downloaded via script).

References
- Human reference proteome (UniProt Proteome ID: `UP000005640`, TaxID 9606), FASTA: `UP000005640_9606.fasta.gz`.
- Gencode annotation: `gencode.v43.basic.annotation.gtf.gz`.

Metabolomics
- Metabolomics Workbench studies used for CSF reference and redox examples:
  - `ST000745` (CSF) — zip and results text
  - `ST001104` — zip and results text
  - `ST002921` — raw zip and results text
- Direct URLs are consolidated in `resources/manifests/metabolomics_urls.txt`.

Additional study supplements
- PXD032767 clinical and run mapping table: Supplementary Table S1
  - Link: https://pmc.ncbi.nlm.nih.gov/articles/instance/9639356/bin/vdac161_suppl_supplementary_table_s1.xlsx
- PXD018301 (Lyden lab) EV matrix: `Human512Reports.xlsx` (place under `data/interim/proteomics/PXD018301/metadata/`). Public direct link varies by source; retrieve from the study’s repository if needed.

One‑Click Retrieval Coverage
- The script `scripts/sh/download_all_data.sh` automates retrieval for:
  - PRIDE datasets (via API): PXD011796, PXD005719, PXD046330, PXD051579, PXD032767, PXD018301
  - Metabolomics Workbench: all entries in `resources/manifests/metabolomics_urls.txt`
  - Xena MET500: expression + clinical
  - GEO GSE12237/GSE125989 series matrices and GPLs
  - GEO GSE186344 supplementary files
  - UniProt FASTA + Gencode GTF

Notes and Limits
- Some proteomics cohorts (e.g., MSV000089062) are used via derived feature tables; full RAW/ID retrieval is not automated here.
- For GSE96860 bigWigs, the volume is large; the script prints instructions to optionally fetch them.
- If a download fails due to network or mirror hiccups, re‑run the script; it uses retry/resume where available.

