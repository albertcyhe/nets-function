Proteomic + Metabolomic DeltaFM Extension Roadmap
===========================================

Context & linkage to core analysis
----------------------------------
- This document extends the existing RNA-centric DeltaFM pipeline (see `README.md`) into proteomic and metabolomic layers.
- Statistical standards mirror the transcriptomic work: paired Wilcoxon with Cliff's delta for within-patient contrasts, Mann-Whitney or linear models for unpaired comparisons, Spearman rho with BH FDR control (alpha=0.05) for correlations, and partial correlations/mixed models controlling the same covariates (neutrophil score, tumor purity, cohort/batch).
- Heavy computations (spectral reanalysis, ssGSEA over large matrices) should run in Colab notebooks stored under `colab/` and synced to Drive, as done previously.

Data inventory & download commands
----------------------------------
Use the reusable PRIDE fetcher (`scripts/sh/fetch_pride_dataset.sh`) to materialise raw/processed proteomic files via `aria2c -c -x 16 -s 16 --max-tries=0 --retry-wait=30`. The script queries the PRIDE v2 API, filters by file category, and emits robust downloads.

| Aim | Dataset | Modality | Notes | Command |
|-----|---------|----------|-------|---------|
| A1  | PXD011796 | NET reference (label-free) | Healthy vs RA/SLE +- stimuli; provides quantitative NET cargo tables. | `scripts/sh/fetch_pride_dataset.sh PXD011796 --categories RAW,RESULT --dest data/raw/PXD011796` |
| A1  | (2025 DIA compendium)* | NET DIA/MS (multi-stimuli) | Use once accession released; plug accession into same script. | `scripts/sh/fetch_pride_dataset.sh <ACCESSION> --categories RAW,RESULT` |
| A1/A2/A3 | PXD005719 | Brain metastatic membrane proteome | Focus on BBB crossing cargo. | `scripts/sh/fetch_pride_dataset.sh PXD005719 --categories RAW,RESULT --dest data/raw/PXD005719` |
| A1/A2/A3 | PXD046330 | HER2+ cells + brain conditioned media secretome | Compare exposed vs control. | `scripts/sh/fetch_pride_dataset.sh PXD046330 --categories RAW,RESULT --dest data/raw/PXD046330` |
| A1/A2/A3 | PXD051579 | HER2+ brain metastasis BBB dysfunction | Time course proteome; contains THBS1 peptides. | `scripts/sh/fetch_pride_dataset.sh PXD051579 --categories RAW,RESULT --dest data/raw/PXD051579` |
| B1/B2 | CSF proteomes (various accessions)** | CSF or plasma proteomics | Use script per accession; prefer datasets with peptide-level evidence. | `scripts/sh/fetch_pride_dataset.sh <CSF_ACCESSION> --categories RESULT --match 'CSF' --dest data/raw/<CSF_ACCESSION>` |
| B1 | CSF metabolomics (Yoo 2017, Ballester 2018, Jang 2025, etc.) | Targeted/untargeted metabolite tables | Direct `studydownload/` endpoints mirrored in `resources/manifests/metabolomics_urls.txt`. | `aria2c -c -x 16 -s 16 --max-tries=0 --retry-wait=30 -d data/raw/metabolomics -i resources/manifests/metabolomics_urls.txt` |

\* Placeholder: update once the 2025 DIA-MS NET resource publishes its PXD accession.

\** For CSF proteomes hosted outside PRIDE (e.g. supplemental XLSX), add URLs to `resources/manifests/csf_proteome_urls.txt` and feed into `aria2c` as shown.

Manifest preparation hints
- For metabolomics tables that lack a manifest, create `resources/manifests/*.txt` with one URL per line. The current manifest lists the `studydownload/` ZIP plus key result tables for ST000745, ST001104, and ST002921.
- To pre-filter PRIDE downloads by filename, pass `--match` (regex). Example: `--match '\\.(raw|mzML)$'` to fetch vendor RAW or mzML files only.
- Run downloads from the project root: `bash scripts/sh/fetch_pride_dataset.sh PXD011796 --categories RESULT`.

Bulk acquisition status (2025-09-17)
------------------------------------
- Colab notebook `00_Bulk_Data_Acquisition.ipynb` now mounts Drive, generates manifests, deduplicates via `SKIP_EXISTING`, and moves files into `/data/raw/<modality>/<accession>/`.
- Transcriptomics: GSE125989, GSE14017, GSE14018, and GSE43837 are stored under `data/raw/transcriptomics/`; GSE184869 provides the curated Excel supplement (series matrix unavailable on GEO).
- Proteomics: PRIDE datasets PXD011796, PXD046330, PXD051579, and PXD005719 were downloaded in full (RAW plus mzID) and reside under `data/raw/proteomics/`.
- Metabolomics: ST000745.zip, ST001104.zip, ST002921_NTM_Rawdata.zip, and their result tables are retrieved via the updated manifest and saved in `data/raw/metabolomics/`.
- Reference: UniProt human reference FASTA (`UP000005640_9606.fasta.gz`) is available under `data/raw/reference/`.

Proteomic preprocessing (Colab notebooks)
----------------------------------------
1. **Notebook templates** (`colab/proteo_*.ipynb`): base them on previous RNA notebooks but mount Drive folders containing PRIDE downloads. Prepared notebooks:
   - `colab/proteo_net_reference.ipynb` – installs msconvert/MSFragger/Philosopher and runs the full pipeline (RAW -> mzML -> MSFragger -> Philosopher) for PXD011796, exporting log2 protein matrices and metadata.
   - `colab/proteo_brain_metastasis.ipynb` – batches the same workflow across PXD005719, PXD046330, and PXD051579, emitting per-dataset `protein_abundance.tsv`/`metadata.tsv` pairs.
   - `colab/proteo_csf.ipynb` (to be drafted) for CSF cohorts.
2. **Dependencies**: notebooks bootstrap OpenJDK and a shared micromamba environment (`pwiz`) that installs ProteoWizard `msconvert` and MSFragger (via bioconda). The install cell provisions the shim `/usr/local/bin/msconvert_pwiz` and automatically resolves the latest `MSFragger-*.jar` from the conda prefix. Add `pyopenms`, `msproteomicstools`, `msqrob2` (via `rpy2`), and `pymzml` for downstream QC as needed.
3. **Raw to quantified matrix**:
   - Convert vendor RAW to mzML using `msconvert` (ProteoWizard CLI). Upload the `.zip` of RAW files to Drive; run `!msconvert *.raw --mzML --filter "peakPicking true 1-"` in Colab.
   - Search strategy:
     - **Conventional quantification (A1/A3)**: run DIA-NN (for DIA compendium) or FragPipe/MaxQuant (for DDA). Export protein groups with LFQ/iBAQ or spectral counts.
     - **Semi-tryptic / N-terminomics (A2)**: MSFragger in semi-tryptic mode (`--search_enzyme_name nonspecific`, `--clip_nTerm=1`) and Philosopher `peptideprophet --nonparam --ppm --expectscore`. Export `.pepXML`/`peptide.tsv` for THBS1 cleavage mapping.
   - Normalize intensities per dataset (median/quantile). Keep log2 scale for downstream Delta metrics.
4. **Annotation harmonisation**: map UniProt accessions to HGNC gene symbols using UniProt release matched to search (store mapping in `resources/annotation/uniprot_<release>.tsv`). For multi-mapping proteins, keep one representative symbol (highest intensity) and store mapping table under `results/tables/<dataset>_protein_gene_map.tsv`.
5. **Batch metadata**: record experimental design (stimulus, disease status, brain vs lung, time point) in `data/processed/<dataset>/<dataset>.metadata.tsv`.
6. **Outputs**:
   - Expression matrix: `data/processed/<dataset>/<dataset>.protein_abundance.tsv` (rows = proteins/genes, columns = samples).
   - Peptide evidence (for degradomics): `data/processed/<dataset>/<dataset>.peptide_level.tsv` plus `thbs1_peptide_events.tsv` enumerating cleavage products.

Proteo-DeltaFM scoring (Aim A1)
---------------------------
1. **Module definition**
   - NET-M (stable marker cargo): start from PXD011796 + DIA compendium. Identify peptides linked to cit-H3/H4/H2A, MPO, PADI4 complexes. Rank proteins by consistency across stimuli (coefficient of variation <0.3 across replicates).
   - NET-F (functional cargo): ELANE, PRTN3, CTSG, NE-associated effectors; include proteins enriched in acute NET release but depleted in chronic/brain metastasis contexts.
   - Store final gene lists in `resources/modules/net_m_proteome.tsv` and `resources/modules/net_f_proteome.tsv` (one gene per line).
2. **Score computation** (run in R `scripts/r/11_proteo_scores.R` or Colab via `GSVA`):
   - For each proteomic dataset, z-standardise log2 intensities per sample.
   - Compute ssGSEA (GSVA param mode) and singscore for both NET-M and NET-F modules.
   - Define `Proteo_DeltaFM = z(score_NET_F) - z(score_NET_M)` (z across samples).
   - Record QC metrics: module coverage (% proteins detected), Cronbach's alpha, pairwise correlations.
3. **Comparisons**
   - Paired contrasts (e.g. conditioned vs control) -> paired Wilcoxon + Cliff's delta; report in `results/tables/<dataset>_proteo_deltafm_paired.tsv`.
   - Brain vs lung metastasis groups -> Mann-Whitney / linear models with neutrophil abundance covariate (proxy from CD66b/MPO intensities).
   - Correlate Proteo-DeltaFM with transcriptomic DeltaFM for overlapping samples (Spearman, partial correlation controlling neutrophil score + purity proxies).
4. **Outputs**: `results/tables/proteo_deltafm_summary.tsv`, `results/figures/F_proteo_deltafm_boxplots.pdf`, `results/figures/F_proteo_deltafm_vs_transcriptome.pdf`.

THBS1 degradomics (Aim A2)
--------------------------
1. **Peptide-centric search**: reuse MSFragger semi-tryptic outputs. Filter peptides mapped to THBS1 with non-tryptic N-termini. Annotate cleavage site position using UniProt isoform-specific sequences.
2. **Cleavage index definition**:
   - `Cleavage_Index = log2(sum intensities of neo-N termini peptides) - log2(sum intensities of full-length tryptic peptides)` per sample (impose minimum peptide count >=2 for each sum).
   - Store results in `data/processed/<dataset>/<dataset>.thbs1_cleavage.tsv`.
3. **Validation**:
   - Cross-check peptides against literature cleavage sites (HtrA1, ADAMTS). Flag novel events.
   - Compare cleavage index between brain and lung cohorts (Mann-Whitney) and correlate with Proteo-DeltaFM / transcriptomic DeltaFM (Spearman, FDR BH).
   - Visualise spectra of top neo-peptides via `pymzml` (save to `results/figures/thbs1_spectrum_<peptide>.pdf`).

Serpin occupancy vs Proteo-DeltaFM (Aim A3)
---------------------------------------
1. **Serpin module**: compile SERPINA3, SERPINB2, SERPINI1, SERPING1, SERPINB1, SERPINB5 etc. Store in `resources/modules/serpin_brain.tsv` (mark astrocyte-enriched serpins with metadata column).
2. **Scoring**: compute ssGSEA/singscore on proteomic matrices to obtain `Serpin_score`.
3. **Analyses**:
   - Correlate `Serpin_score` with `Proteo_DeltaFM` (expect negative correlation in brain metastasis proteomes). Report rho, p, FDR, 95% CI in `results/tables/proteo_serpin_assoc.tsv`.
   - Partial correlations controlling neutrophil proxies and astrocyte markers (GFAP, S100B) to isolate effect.
   - Overlay on transcriptomic Serpin scores (GSE184869) for matched subjects.

Metabolomic extensions (B1 & B2)
----------------------------------
1. **Data ingestion**: `colab/metabolomics_preprocessing.ipynb` unpacks Workbench `studydownload` archives (`ST000745.zip`, `ST001104.zip`, `ST002921_NTM_Rawdata.zip`) into `data/processed/metabolomics/raw_tables/`, harmonising HMDB IDs when provided.
2. **Redox metrics (B1)**:
   - Notebook computes GSH/GSSG, NADPH/NADP+, MetSO/Met, and lactate/pyruvate ratios per sample; output saved to `data/processed/metabolomics/redox_ratios.tsv`.
   - Correlate each ratio with transcriptomic DeltaFM and Proteo-DeltaFM within matching cohorts (`Spearman`, BH FDR). Store in `results/tables/metabolomics_redox_corr.tsv`.
   - Group comparisons (brain vs non-brain): Mann-Whitney + Cliff's delta.
3. **MPO chlorination traces (B2)**:
   - In proteomics: search for 3-chlorotyrosine modifications (add variable mod in MSFragger: `Chlorination on Y, mass shift +34.9694 Da`). Summarise detection frequency and intensities; define `MPO_trace_score` per sample.
   - In metabolomics: identify chloramine-related compounds via HMDB cross-reference; compute abundance z-scores.
   - Compare `MPO_trace_score` across organ sites and correlate with Proteo-DeltaFM / Serpin scores.

Cross-layer integration & modelling
-----------------------------------
1. **Sample alignment**: build `results/tables/multilayer_sample_map.tsv` linking samples across transcriptome, proteome, degradome, metabolome (columns: cohort, organ, patient ID, sample type, available layers).
2. **Derived metrics table**: merge DeltaFM (RNA), Proteo-DeltaFM, Serpin score, THBS1 cleavage index, ROS metrics, MPO trace into `results/tables/multilayer_features.tsv`.
3. **Mixed-effects model** (`scripts/r/12_mixed_effects.R`):
   - Formula: `Functional_proxy ~ DeltaFM * Organ + Proteo_DeltaFM + Serpin_score + ROS_ratio + (1 | Cohort)`.
   - Fit using `lme4::lmer` with Satterthwaite df (via `lmerTest`). Report fixed effects, interaction (DeltaFM x organ), conditional R^2.
   - Evaluate organ-specific slopes via emmeans (`emmeans::emtrends`).
   - Diagnostic plots: residual vs fitted, Q-Q normal plot, influence (Cook's distance) saved to `results/figures/mixed_effects_diagnostics.pdf`.
4. **Visual storytelling**: construct the "double sandwich" schematic summarising Proteo-DeltaFM and THBS1 cleavage divergences between brain vs lung metastasis (`results/figures/F_double_sandwich.pdf`).

Directory conventions & deliverables
------------------------------------
- **Raw downloads**: `data/raw/<dataset>/` (tar/RAW/mzML/metadata). Keep `checksums.sha256` after download.
- **Processed matrices**: `data/processed/<dataset>/` (protein and metabolite tables, cleavage indices, covariates).
- **Resources**: `resources/modules/` (module gene lists), `resources/manifests/` (aria2 manifests), `resources/annotation/` (UniProt <-> gene symbol maps).
- **Scripts**:
  - `scripts/sh/fetch_pride_dataset.sh` - PRIDE downloader.
  - Add Colab notebooks under `colab/` with matching snake_case names; document executed notebook per run in `logs/colab_runs.md`.
  - Future R scripts: `scripts/r/11_proteo_scores.R`, `scripts/r/12_mixed_effects.R`.
- **Results**: `results/tables/` and `results/figures/` for each aim (use prefixes `proteo_`, `metabolomics_`, `mixed_`).

Next actions checklist
----------------------
- [ ] Populate `resources/modules/net_m_proteome.tsv`, `net_f_proteome.tsv`, and `serpin_brain.tsv` from PXD011796 + literature.
- [x] Draft Colab notebook `proteo_net_reference` with RAW->mzML->MSFragger workflow; export protein-level table.
- [x] Draft Colab notebook `proteo_brain_metastasis` covering PXD005719/PXD046330/PXD051579.
- [ ] Run `scripts/sh/fetch_pride_dataset.sh PXD046330 --categories RESULT` to obtain conditioned media quant tables; compute Proteo-DeltaFM.
- [ ] Configure MSFragger semi-tryptic searches for THBS1 cleavage discovery (store parameter file in `resources/msfragger/semi_tryptic.params`).
- [x] Assemble metabolomics manifests (Yoo 2017 etc.).
- [x] Draft Colab notebook `metabolomics_preprocessing` for archive extraction + ratio computation.
- [ ] Compute redox ratios once processed metabolite matrices are harmonised.
- [ ] Integrate layers into mixed-effects model and craft publishable figures.
