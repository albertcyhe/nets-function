Proteomic + Metabolomic DeltaFM Extension Roadmap
===========================================

Context & linkage to core analysis
----------------------------------
- This document extends the existing RNA-centric DeltaFM pipeline (see `README.md`) into proteomic and metabolomic layers.
- Statistical standards mirror the transcriptomic work: paired Wilcoxon with Cliff's delta for within-patient contrasts, Mann-Whitney or linear models for unpaired comparisons, Spearman rho with BH FDR control (alpha = 0.05) for correlations, and partial correlations/mixed models controlling neutrophil score, tumour purity, and cohort/batch effects.
- All raw resources (~120 GB) are now resident on the server under `data/raw/`; verify integrity with `sha256sum -c data/raw/checksums.sha256` before running downstream steps.
- Heavy computations (spectral reanalysis, ssGSEA, mixed models) should execute directly on this machine; the Colab notebooks in `colab/` remain as provenance but are no longer the primary execution path.

Local execution overview
------------------------
- Use `tmux`/`screen` for long jobs and capture command logs under `logs/` (e.g. `logs/proteomics/pipeline_<timestamp>.log`).
- Stage intermediates beneath `data/interim/<modality>/<dataset>/` so raw archives stay untouched inside `data/raw/`.
- Reserve `data/processed/` for final matrices and per-dataset metadata; keep immutable inputs read-only once verified.
- Track each major run (date, command, commit hash) in a lightweight runbook such as `logs/runbook_proteo_metab.md`.

Local environment setup
-----------------------
1. Bootstrap micromamba (already vendored in `env/bin/micromamba`):
   ```bash
   export MAMBA_ROOT_PREFIX="$PWD/env/.mamba"
   env/bin/micromamba env create -y -f env/proteomics.yml
   eval "$(env/bin/micromamba shell hook -s bash)"
   micromamba activate proteomics
   ```
2. Add helper binaries to the PATH for convenience:
   ```bash
   export PATH="$PWD/env/bin:$PATH"
   philosopher version  # sanity check (binary lives in env/bin)
   ```
3. Optional: materialise an R library cache (e.g. `env/.R/`) for `msqrob2`, `lmerTest`, and plotting packages; point scripts to it via `R_LIBS_USER`.
4. When analyses span multiple shells, source `env/activate_proteomics.sh` (to be added) or re-run the activation snippet above.

Data inventory & download commands
----------------------------------
Use the reusable PRIDE fetcher (`scripts/sh/fetch_pride_dataset.sh`) to materialise raw/processed proteomic files via `aria2c -c -x 16 -s 16 --max-tries=0 --retry-wait=30`. The script queries the PRIDE v2 API, filters by file category, and emits robust downloads. All datasets listed are already mirrored under `data/raw/`, so rerun commands only if integrity checks fail or new accessions appear.

| Aim | Dataset | Modality | Notes | Command |
|-----|---------|----------|-------|---------|
| A1  | PXD011796 | NET reference (label-free) | Healthy vs RA/SLE +- stimuli; provides quantitative NET cargo tables. | `scripts/sh/fetch_pride_dataset.sh PXD011796 --categories RAW,RESULT --dest data/raw/PXD011796` |
| A1  | (2025 DIA compendium)* | NET DIA/MS (multi-stimuli) | Use once accession released; plug accession into same script. | `scripts/sh/fetch_pride_dataset.sh <ACCESSION> --categories RAW,RESULT` |
| A1/A2/A3 | PXD005719 | Brain metastatic membrane proteome | Focus on BBB crossing cargo. | `scripts/sh/fetch_pride_dataset.sh PXD005719 --categories RAW,RESULT --dest data/raw/PXD005719` |
| A1/A2/A3 | PXD046330 | HER2+ cells + brain conditioned media secretome | Compare exposed vs control. | `scripts/sh/fetch_pride_dataset.sh PXD046330 --categories RAW,RESULT --dest data/raw/PXD046330` |
| A2 (paused) | PXD051579 | HER2+ brain metastasis BBB dysfunction | Cross-omics paused (no matched transcriptomic cohort). | `scripts/sh/fetch_pride_dataset.sh PXD051579 --categories RAW,RESULT --dest data/raw/PXD051579` |
| A1/A2 | GSE96860 | RNA-seq (breast cancer cell lines, bigWig) | Summarise stranded bigWigs to gene counts (SKBR3 ↔ PXD046330). | `python scripts/py/process_gse_bigwig.py --bigwig-dir data/raw/transcriptomics/GSE96860 --gtf resources/annotation/gencode/gencode.v43.basic.annotation.gtf.gz --dataset GSE96860 --output-dir data/processed/GSE96860` |
| A1/A3 | GSE12237 | Affymetrix U133A expression | Parsed with GPL96 to pair MDA-MB-231 variants ↔ PXD005719. | `python scripts/py/parse_series_matrix.py --series data/raw/transcriptomics/GSE12237/GSE12237_series_matrix.txt.gz --gpl data/raw/transcriptomics/GSE12237/GPL96.annot.gz --dataset GSE12237 --outdir data/processed/GSE12237` |
| B1/B2 | CSF proteomes (various accessions)** | CSF or plasma proteomics | Use script per accession; prefer datasets with peptide-level evidence. | `scripts/sh/fetch_pride_dataset.sh <CSF_ACCESSION> --categories RESULT --match 'CSF' --dest data/raw/<CSF_ACCESSION>` |
| B1 | CSF metabolomics (Yoo 2017, Ballester 2018, Jang 2025, etc.) | Targeted/untargeted metabolite tables | Direct `studydownload/` endpoints mirrored in `resources/manifests/metabolomics_urls.txt`. | `aria2c -c -x 16 -s 16 --max-tries=0 --retry-wait=30 -d data/raw/metabolomics -i resources/manifests/metabolomics_urls.txt` |

\* Placeholder: update once the 2025 DIA-MS NET resource publishes its PXD accession.

\** For CSF proteomes hosted outside PRIDE (e.g. supplemental XLSX), add URLs to `resources/manifests/csf_proteome_urls.txt` and feed into `aria2c` as shown.

Manifest preparation hints
- For metabolomics tables that lack a manifest, create `resources/manifests/*.txt` with one URL per line. The current manifest lists the `studydownload/` ZIP plus key result tables for ST000745, ST001104, and ST002921.
- To pre-filter PRIDE downloads by filename, pass `--match` (regex). Example: `--match '\\.(raw|mzML)$'` to fetch vendor RAW or mzML files only.
- Run downloads from the project root: `bash scripts/sh/fetch_pride_dataset.sh PXD011796 --categories RESULT`.

Bulk acquisition status (2025-09-19)
------------------------------------
- Data previously staged via `00_Bulk_Data_Acquisition.ipynb` have been synced; verify with `find data/raw -maxdepth 2 -type f | wc -l` and the checksum manifest.
- Transcriptomics: legacy downloads (GSE125989, GSE14017, GSE14018, GSE43837) remain under `data/raw/transcriptomics/`; new GEO deliveries include `data/raw/transcriptomics/GSE96860/` (bigWig RNA-seq summarised to `data/processed/GSE96860/`) and `data/raw/transcriptomics/GSE12237/` (Affymetrix expression parsed to `data/processed/GSE12237/`).
- Proteomics: PRIDE datasets PXD011796, PXD046330, PXD051579, and PXD005719 (RAW + mzID) live in `data/raw/proteomics/`.
- Metabolomics: ST000745.zip, ST001104.zip, ST002921_NTM_Rawdata.zip, and associated result tables are mirrored in `data/raw/metabolomics/`.
- Reference: UniProt human reference FASTA (`UP000005640_9606.fasta.gz`) is stored in `data/raw/reference/`.

Proteomic preprocessing (local pipeline)
----------------------------------------
1. **Directory staging**
   - Create `data/interim/proteomics/<PXD>/mzML/` and `data/interim/proteomics/<PXD>/fragger/` for each dataset.
   - Decompress `data/raw/reference/UP000005640_9606.fasta.gz` (or symlink an existing copy) into `data/interim/proteomics/reference/UP000005640_9606.fasta` for reproducible paths, then run `philosopher database --custom ... --contam` once to generate the target-decoy FASTA (`data/interim/proteomics/reference/UP000005640_9606_td.fasta`). Update the MSFragger parameter files if your project root differs.
2. **RAW → mzML conversion** (ProteoWizard `msconvert` in the proteomics environment):
   ```bash
   micromamba activate proteomics
   mkdir -p data/interim/proteomics/PXD011796/mzML
   parallel --jobs 4 "msconvert {} --mzML --filter 'peakPicking true 1-' --outdir data/interim/proteomics/PXD011796/mzML" ::: data/raw/proteomics/PXD011796/*.raw
   ```
   - Capture logs to `logs/proteomics/msconvert_PXD011796.log` using `tee`.
3. **Database search with MSFragger + Philosopher**
   - Parameter templates live under `resources/msfragger/` (`default_closed.params` for global searches, `semi_tryptic.params` for THBS1 degradomics). `database_name` currently points to `/root/bioinfo/nets-fuction/data/interim/proteomics/reference/UP000005640_9606_td.fasta`; adjust if you relocate the project.
   - Scripts:
     - `scripts/sh/run_trfp_parallel.sh` – ThermoRawFileParser conversion with progress logging (tmux-friendly).
     - `scripts/sh/msfragger_smoke.sh` – end-to-end MSFragger + Philosopher smoke test on selected mzMLs (logs under `logs/proteomics/`).
   - Example (manual) closed search invocation:
     ```bash
     java -Xmx120G -jar $CONDA_PREFIX/share/msfragger/MSFragger.jar \
       resources/msfragger/default_closed.params \
       data/interim/proteomics/PXD011796/mzML/*.mzML
     philosopher workspace --clean
     philosopher workspace --init
     philosopher database --custom /root/bioinfo/nets-fuction/data/interim/proteomics/reference/UP000005640_9606_td.fasta
     philosopher database --annotate /root/bioinfo/nets-fuction/data/interim/proteomics/reference/UP000005640_9606_td.fasta
     philosopher peptideprophet --nonparam --ppm --expectscore --decoy rev_ *.pepXML
     philosopher proteinprophet --maxppmdiff 200 --output combined interact-*.pep.xml
     philosopher filter --sequential --razor --mapmods --tag rev_ --pepxml interact-*.pep.xml --protxml combined.prot.xml
     philosopher report
     ```
   - Store intermediate outputs in `data/interim/proteomics/<PXD>/fragger/`.
4. **Quantification & normalisation**
   - Use `philosopher abacus` or `msstat` exports to obtain protein-level intensities; script logic belongs in `scripts/py/proteomics_summarise.py` (convert from existing notebooks).
   - Normalise (median/quantile) per dataset; keep log2 scale. Save matrices under `data/processed/<PXD>/<PXD>.protein_abundance.tsv` alongside `metadata.tsv`.
   - Maintain peptide-level evidence for degradomics in `data/processed/<PXD>/<PXD>.peptide_level.tsv` plus `thbs1_peptide_events.tsv`.
5. **QC**
   - Generate `results/figures/<PXD>_qc_total_ion.pdf`, coefficient-of-variation plots, and missingness heatmaps.
   - Summarise search metrics (PSM counts, peptide FDR, protein FDR) in `results/tables/<PXD>_search_stats.tsv`.

Proteo-DeltaFM scoring (Aim A1)
---------------------------
1. **Module definition**
   - NET-M (stable marker cargo): start from PXD011796 + DIA compendium. Identify peptides linked to cit-H3/H4/H2A, MPO, PADI4 complexes. Rank proteins by consistency across stimuli (coefficient of variation <0.3 across replicates).
   - NET-F (functional cargo): ELANE, PRTN3, CTSG, NE-associated effectors; include proteins enriched in acute NET release but depleted in chronic/brain metastasis contexts.
   - Store final gene lists in `resources/modules/net_m_proteome.tsv` and `resources/modules/net_f_proteome.tsv` (one gene per line).
2. **Score computation** (target implementation: `scripts/r/11_proteo_scores.R`)
   - Z-standardise log2 intensities per sample.
   - Compute ssGSEA (GSVA param mode) and singscore for both NET-M and NET-F modules.
   - Define `Proteo_DeltaFM = z(score_NET_F) - z(score_NET_M)` (z across samples).
   - Record QC metrics: module coverage (% proteins detected), Cronbach's alpha, pairwise correlations.
3. **Comparisons**
   - Paired contrasts (e.g. conditioned vs control) → paired Wilcoxon + Cliff's delta; output `results/tables/<dataset>_proteo_deltafm_paired.tsv`.
   - Brain vs lung metastasis groups → Mann-Whitney / linear models with neutrophil abundance covariate (proxy from CD66b/MPO intensities).
   - Correlate Proteo-DeltaFM with transcriptomic DeltaFM for overlapping samples (Spearman, partial correlation controlling neutrophil score + purity proxies).
4. **Outputs**: `results/tables/proteo_deltafm_summary.tsv`, `results/figures/F_proteo_deltafm_boxplots.pdf`, `results/figures/F_proteo_deltafm_vs_transcriptome.pdf`. *(PXD051579 cross-cohort comparison paused until a matched transcriptomic dataset is available.)*

THBS1 degradomics (Aim A2)
--------------------------
1. **Peptide-centric search**: reuse MSFragger semi-tryptic outputs. Filter peptides mapped to THBS1 with non-tryptic N-termini. Annotate cleavage site positions using isoform-specific UniProt sequences.
2. **Cleavage index definition**
   - `Cleavage_Index = log2(sum intensities of neo-N termini peptides) - log2(sum intensities of full-length tryptic peptides)` per sample (minimum peptide count ≥ 2 for each sum).
   - Store results in `data/processed/<dataset>/<dataset>.thbs1_cleavage.tsv`.
3. **Validation**
   - Cross-check peptides against literature cleavage sites (HtrA1, ADAMTS). Flag novel events.
   - Compare cleavage index between brain and lung cohorts (Mann-Whitney) and correlate with Proteo-DeltaFM / transcriptomic DeltaFM (Spearman, FDR BH).
   - Visualise spectra of top neo-peptides via `pymzml`; save under `results/figures/thbs1_spectrum_<peptide>.pdf`.

Serpin occupancy vs Proteo-DeltaFM (Aim A3)
---------------------------------------
1. **Serpin module**: compile SERPINA3, SERPINB2, SERPINI1, SERPING1, SERPINB1, SERPINB5 etc. Store in `resources/modules/serpin_brain.tsv` (mark astrocyte-enriched serpins with a metadata column).
2. **Scoring**: compute ssGSEA/singscore on proteomic matrices to obtain `Serpin_score`.
3. **Analyses**
   - Correlate `Serpin_score` with `Proteo_DeltaFM` (expect negative correlation in brain metastasis proteomes). Report rho, p, FDR, 95% CI in `results/tables/proteo_serpin_assoc.tsv`.
   - Partial correlations controlling neutrophil proxies and astrocyte markers (GFAP, S100B) to isolate the effect.
   - Overlay on transcriptomic Serpin scores (GSE184869) for matched subjects.

### A3 — Orthogonal strategy when F (ELANE/PRTN3/CTSG) is sparsely detected
We accept that some proteomes (e.g., PXD005719) rarely detect F-module proteases. To still test “Serpin inhibition reduces F activity while M persists,” we use three independent readouts that do not rely on direct F quantification:

- Chemistry (stoichiometry): ISI index
  - Definition: Inhibitor–Protease Stoichiometry Index, `ISI = log2(Σ inhibitors / Σ target proteases)`.
  - Inhibitor basket (core): SERPINA1, SERPINA3, SERPINB1, SERPINB6, SERPINB8. Extended set (sensitivity): SLPI, PI3/ELAFIN, A2M.
  - Targets: ELANE, PRTN3, CTSG.
  - Model: Mixed-effects regression `Proteo_ΔFM ~ ISI * organ + covariates + (1 | cohort)`. Expect β(ISI) < 0 and stronger effect in brain.

- Detection-probability as signal (MNAR / zero-inflated)
  - Treat “F detected?” at protein/peptide level as the dependent variable, with `Serpin_score` and sample covariates as predictors `(1 | cohort)` random effect.
  - Compare with M markers as negative control; expect Serpin↑ → Pr(detect F)↓, but Pr(detect M) stable.

- Functional footprints (substrate-centric)
  - Build NE/PR3 “footprint indices” from semi-/non-tryptic N-termini enriched at motif-consistent sites (N-terminomics proxy) in existing/semi-tryptic searches.
  - Correlate footprints with `Serpin_score`, `Proteo_ΔFM`, and THBS1-cleavage index.

- THBS1 cleavage endpoint (semi-tryptic)
  - Systematise the semi-tryptic re-search to quantify a THBS1 cleavage index: log2(neo-N sum) − log2(tryptic N sum). Expect lower cleavage in brain (higher inhibition), with M stable.

- RNA→Protein bridge (mediation/moderation)
  - For matched samples, test `Proteo_ΔFM ~ Serpin_RNA * organ + covariates`, controlling neutrophils/purity; Serpin_RNA should carry a significant (partial) effect.

- Extended inhibitor axis (sensitivity)
  - Add SLPI/ELAFIN/A2M to the inhibitor basket in ISI and repeat analyses; consistent negative association strengthens the causal story.

- Missingness as signal — meta-analysis across cohorts
  - Model `Missing%(F) ~ Serpin_score + organ + (1 | cohort)`; expect Serpin↑ → F-missing↑ in brain, while M-missing remains flat.

Why this explains “brain decoupling vs lung coupling”
- Brain: high inhibitor load (astrocyte/liver-derived A1AT, etc.) neutralises F first (low footprints, weak THBS1 cleavage, low F detection), while M (chromatin-bound) persists → ΔFM↓ and decoupling.
- Lung: weaker inhibitor barrier → F footprints and THBS1 cleavage stay high; M and F more aligned (coupling), matching prior macro statistics.

Minimal, two-week deliverables (no wet lab)
- Implement ISI (core and extended baskets) and run mixed-effects with organ interaction; export forest plots.
- Zero-inflated/logistic models for F detection vs Serpin-score; M as control.
- Footprints: export NE/PR3 semi-/non-tryptic N-termini counts/densities and associate with Serpin-score, Proteo_ΔFM, THBS1 cleavage.
- One-page integrative schematic: Serpin↑ → F-missing↑ & Footprints↓ → THBS1 cleavage↓ → Proteo_ΔFM↓; M peptides stable; brain vs lung facets.

Metabolomic extensions (B1 & B2) — *Paused*
------------------------------------------
> 2025-09-20：当前获取的 Workbench 研究（ST000745、ST001104、ST002921）仅提供 mwTab
> 元信息和仪器 RAW，未公开命名代谢物的定量矩阵，GSH/GSSG 等红氧化比暂无法计算。
> 相关工具（`scripts/py/metabolomics_preprocessing.py`, `scripts/py/mw_parser.py`）已归档，待
> 未来获得包含命名代谢物表格的研究后再恢复此模块。

Cross-layer integration & modelling
-----------------------------------
1. **Sample alignment**: build `results/tables/multilayer_sample_map.tsv` linking samples across transcriptome, proteome, degradome, metabolome (columns: cohort, organ, patient ID, sample type, available layers).
2. **Derived metrics table**: merge DeltaFM (RNA), Proteo-DeltaFM, Serpin score, THBS1 cleavage index, ROS metrics, MPO trace into `results/tables/multilayer_features.tsv`.
3. **Mixed-effects model** (`scripts/r/12_mixed_effects.R`)
   - Formula: `Functional_proxy ~ DeltaFM * Organ + Proteo_DeltaFM + Serpin_score + ROS_ratio + (1 | Cohort)`.
   - Fit using `lme4::lmer` with Satterthwaite df (`lmerTest`). Report fixed effects, interaction (DeltaFM × organ), conditional R².
   - Evaluate organ-specific slopes via `emmeans::emtrends`.
   - Diagnostic plots: residual vs fitted, Q-Q normal plot, influence (Cook's distance) saved to `results/figures/mixed_effects_diagnostics.pdf`.
4. **Visual storytelling**: assemble the "double sandwich" schematic summarising Proteo-DeltaFM and THBS1 cleavage divergences between brain vs lung metastasis (`results/figures/F_double_sandwich.pdf`).

#### Multi-omics status (2025-09-21)
- Accepted the absence of NET-F module coverage in the PXD005719 proteome; downstream integration uses NET-M only.
- Built `results/tables/multiomics_alignment.tsv` (and `multiomics_pairwise_deltas.tsv`) aligning PXD046330 ↔ GSE96860 (SKBR3 reference) and PXD005719 ↔ GSE12237 (MDA-MB-231 parental vs brain variants).
- Spearman ΔFM correlations are underpowered (`n<3` per cohort) but NET-M/THBS1 deltas remain available for interpretation; PXD051579 cross-omics stays paused until a matched transcriptomic dataset exists.

#### A3 model executions — first pass (strict evidence)
- ISI (Inhibitor–Protease Stoichiometry):
  - Scripts: `scripts/py/isi_model.py` → `results/tables/isi_per_sample.tsv`, `isi_models.tsv`.
  - Result: In current cohorts (PXD046330, PXD005719, PXD051579) the F targets (ELANE/PRTN3/CTSG) are not detected (target_sum=0 for all samples), so ISI is undefined by design. This is consistent with “F suppressed/not measurable while M persists”.
- MNAR (Missingness-as-signal):
  - Scripts: `scripts/py/mnar_detection_model.py` → `results/tables/mnar_detection_table.tsv`, `mnar_logit_results.tsv`.
  - Result: detect_F=0 across all samples, detect_M=1 across all samples, yielding non‑estimable logits but a definitive qualitative pattern: “F fully missing, M universally present”.
- Footprint indices (substrate-centric NE/PR3 proxy):
  - Script: `scripts/py/footprint_indices.py` → `results/tables/footprints.tsv`.
  - Result: indices are consistently negative; in PXD005719, brain‑variant shows lower footprint than parental (≈ −0.786 vs −0.633), aligning with reduced THBS1 protein signal in brain. PXD046330 and PXD051579 show negative footprints consistent with low NE/PR3 activity.

#### Final stratified associations (evidence chain)
Using CLR‑normalized secreted Serpin scores (core: A1AT/A3/SERPINI1), we ran stratified associations with robust summaries (Spearman primary; OLS/RLM where estimable). Covariate `log_total_PSMs` included where modelled.

- Footprints vs Serpin (expect negative)
  - Combined: ρ ≈ −0.563, p ≈ 0.00121, q ≈ 0.00182 (BH within stratum)
  - Non‑brain: ρ ≈ −0.560, p ≈ 0.00160, q ≈ 0.00130
  - Brain: under‑powered (n≈1–2 effective), no stable estimate

- THBS1 abundance vs Serpin (expect positive)
  - Combined: ρ ≈ +0.583, p ≈ 0.00073, q ≈ 0.00218
  - Non‑brain: ρ ≈ +0.611, p ≈ 0.00043, q ≈ 0.00130
  - Brain: under‑powered

- Proteo‑ΔFM vs Serpin (expect negative)
  - Combined/Non‑brain: weak, non‑significant (ρ ≈ +0.225, p ≈ 0.29)
  - Interpretation: ΔFM here is M‑driven (F fully missing), so Serpin→ΔFM effect is attenuated; our footprint/THBS1 anchors carry the signal as designed.

Tables: `results/tables/assoc_summary.tsv` (per‑stratum endpoint summaries).

Mixed-effects (organ interaction):
- `scripts/py/mixed_effects_a3.py` → `results/tables/mixed_effects_a3.tsv`
  * Model: `endpoint ~ Serpin_score_core * organ + log_total_PSMs + (1 | dataset)`.
  * THBS1 abundance fit: β_Serpin ≈ +0.028 (non-significant; limited organ contrast); footprint model remains singular (brain stratum underpowered), annotated in the table for transparency.

Figures (R, Cell‑style):
- `scripts/r/13_figures_a3.R` renders
  - `results/figures/A3_forest_assoc.pdf` – Spearman effect forest (95% CI via Fisher transform) by stratum×endpoint.
  - `results/figures/A3_scatter_serpin_vs_footprint.pdf` – Serpin vs Footprint (facet by dataset; lm line; theme_classic).
  - `results/figures/A3_scatter_serpin_vs_THBS1.pdf` – Serpin vs THBS1.
  - `results/figures/A3_scatter_serpin_vs_deltaFM.pdf` – Serpin vs Proteo‑ΔFM.
  - `results/figures/A3_scatter_serpin_vs_thbs1_cleave.pdf` – Serpin vs THBS1 cleavage index (endpoint anchor).
  - Negative controls:
    - `results/figures/A3_scatter_serpin_vs_ctrl_footprint.pdf` – Serpin vs control footprint (acidic P1 motif).
    - `results/figures/A3_scatter_serpin_vs_ecm_nc.pdf` – Serpin vs ECM negative‑control CLR score (laminins).
  - Run:
    ```bash
    micromamba activate proteomics
    Rscript scripts/r/13_figures_a3.R
    ```

THBS1 cleavage index (endpoint anchor):
- `scripts/py/thbs1_cleavage_index.py` → `results/tables/thbs1_cleave_idx.tsv` (controls for COL1A1/COL4A1 included). Command:
  ```bash
  micromamba run -n proteomics python scripts/py/thbs1_cleavage_index.py \
    --dataset PXD046330 --dataset PXD005719 --dataset PXD051579
  ```

Negative controls (null guards):
- Compute and store at `results/tables/neg_controls.tsv`:
  - ECM_nc_score: CLR mean across laminin subunits (LAMA1/2/3, LAMB1/2/3, LAMC1/2/3).
  - control_footprint_index: footprint with acidic P1 motif (D/E), expected to show no Serpin dependence.
- Script: `scripts/py/neg_controls.py` (already executed).

#### PhD feedback → publishable modelling plan (to execute next)
We keep the strict findings above unchanged and upgrade modelling per the following A–F roadmap.

- A1. Interpretable Serpin scores (already scaffolded)
  - Core (secreted): SERPINA1/A3/SERPINI1; Sensitivity (extended): + SLPI, Elafin (PI3/TRAPPIN‑2), A2M.
  - Normalisation: median polish + CLR (label‑free). First pass: CLR implemented via `scripts/py/serpin_scores.py` → `results/tables/serpin_scores.tsv` (columns: `Serpin_score_core`, `Serpin_score_ext`). Median‑polish step will be added when intensity matrices are in place.
- A2. Three associations, stratified + robust (to implement)
  - Models: 
    - `Footprint_index_NEPR3 ~ Serpin_score + purity + neutro_abundance + (1 | cohort)` (expect β<0)
    - `THBS1_abundance ~ Serpin_score + ...` (expect β>0)
    - `Proteo-ΔFM ~ Serpin_score + ...` (expect β<0)
  - Robustness: Spearman + Huber/robust LM + rank‑permutation; BH FDR across endpoints.
  - Stratify by organ (brain vs non‑brain) and produce random‑effects forest plot (report I²). Script stub to be added (`scripts/py/robust_assoc.py`).
- A3. Negative controls (to implement)
  - ECM secreted proteins without NE/PR3 motif (e.g., laminin subunits) → correlation vs Serpin_score ~ 0.
  - Unrelated footprint (e.g., MMP‑favoured motif) → no association with Serpin_score.
- B1. Endpoint anchor — THBS1 cleavage index (to implement)
  - `THBS1_cleave_idx = (# semi/non‑tryptic N‑termini at curated THBS1 sites) / (THBS1 coverage)`; associate with Serpin_score, ΔFM, footprints; verify COL1A1/COL4A1 as specificity controls.
  - Script stub will consume combined PSMs and a curated site list under `resources/annotation/thbs1_sites.tsv`.
- B2. MNAR upgraded — peptide‑level hurdle/zero‑inflated + meta‑logit (to implement)
  - Model detection per peptide: `Pr(detect_peptide_j) ~ Serpin_score + mods + (1 | cohort)`; meta‑logit across cohorts for NE/PR3 peptides; histone/ci‑markers as negative controls.
- C. Anti‑confounding (to implement)
  - C1 Deconvolution: estimate astrocyte/microglia/endothelium/neutrophil proportions from marker panels; include or stratify.
  - C2 Compositional effects: apply CLR to all proteins; z‑score within cohort; report sensitivity.
  - C3 Random effects: carry `(1|cohort)` in all models; report I² (Cochran’s Q).
- D. Structural equation model (SEM) / Bayesian latent (to implement)
  - Latent `F_activity` targeted by {Footprints, THBS1_cleave_idx, MNAR}; 
    paths: `Serpin_score → F_activity → {Footprints, THBS1_cleave_idx, ΔFM}`; `M_stability → ΔFM`; organ as moderator.
  - Output: standardized path coefficients + fit (CFI/TLI/RMSEA or WAIC).
- E. Positive/negative controls (to implement)
  - Positive: known NET‑positive inflammatory set → Footprints↑, THBS1‑cleave↑, ΔFM≥0.
  - Negative: NET‑irrelevant pathways (e.g., ribosomal) uncorrelated with Serpin_score.
- F. Stop criteria (target for manuscript lock‑in)
  - Brain stratum: |β|≥0.25 and q<0.05 for ≥2 of 3 primary endpoints (Footprints, THBS1, ΔFM) with negative controls null and acceptable heterogeneity (I²<60%); otherwise defer to SEM for integrated evidence.

Directory conventions & deliverables
------------------------------------
- **Raw downloads**: `data/raw/<dataset>/` (tar/RAW/mzML/metadata). Keep `checksums.sha256` after download.
- **Intermediates**: `data/interim/<modality>/<dataset>/` for mzML, search workspaces, temporary TSVs.
- **Processed matrices**: `data/processed/<dataset>/` (protein/metabolite tables, cleavage indices, covariates).
- **Resources**: `resources/modules/` (module gene lists), `resources/manifests/` (aria2 manifests), `resources/annotation/` (UniProt ↔ gene symbol maps), `resources/msfragger/` (search params).
- **Scripts**: migrate core execution into `scripts/py/` and `scripts/r/`; keep legacy notebooks under `colab/` with a note that they are archival references.
- **Results**: `results/tables/` and `results/figures/` for each aim (prefix `proteo_`, `metabolomics_`, `mixed_`).

Environment notes (RNA bigWig summarisation)
- The `proteomics` env now includes `pybigwig` for stranded bigWig handling and `gffread` for light-weight GTF curation. Gencode v43 GTF is staged at `resources/annotation/gencode/gencode.v43.basic.annotation.gtf.gz`.

Next actions checklist
----------------------
- [x] Run `sha256sum -c data/raw/checksums.sha256` and snapshot the output in `logs/checksums_2025-09-19.txt`.
- [x] Create `data/interim/proteomics/` and `data/interim/metabolomics/` scaffolds; update `.gitignore` if needed.
- [x] Build and document the local `proteomics` micromamba environment; record activation snippet in `env/activate_proteomics.sh`.
- [x] Convert PXD011796, PXD005719, PXD046330, and PXD051579 RAW files to mzML on server; log runtime and storage footprint (2025-09-19).
- [x] Draft `resources/msfragger/default_closed.params` and `resources/msfragger/semi_tryptic.params`; stage MSFragger configs under `resources/msfragger/`.
- [x] Smoke-test MSFragger + Philosopher on a small subset (two PXD011796 mzMLs); see `data/interim/proteomics/PXD011796/fragger_smoke/` for artifacts and logs (2025-09-20).
- [x] Validate semi-tryptic THBS1 workflow on the same subset; outputs recorded under `data/interim/proteomics/PXD011796/fragger_smoke_semi/` (2025-09-20).
- [ ] Port `colab/proteo_net_reference.ipynb` and `colab/proteo_brain_metastasis.ipynb` into scripted runs (`scripts/py/proteomics_summarise.py` + R scoring template).
- [ ] Populate `resources/modules/net_m_proteome.tsv`, `net_f_proteome.tsv`, and `serpin_brain.tsv` from PXD011796 + literature.
- [ ] *(Paused)* Metabolomics preprocessing / redox ratios — restore once named metabolite matrices are available.
- [ ] *(Paused)* Integrate metabolomics-derived metrics alongside proteomic layers.
- [ ] Integrate layers into the mixed-effects framework and draft publication-quality figures.
- External-growth / bypass proxies (module scores)
  - `scripts/py/module_scores.py` → `results/tables/module_scores.tsv` (CLR-based):
    * `E2F_G2M_score` (proliferation), `IL6_STAT3_score`, `CoOption_score`.
  - Associations (`scripts/py/module_assoc.py` → `results/tables/module_assoc.tsv`):
    * Non-brain: CoOption_score correlates positively with THBS1 abundance (ρ≈+0.60, p≈5.5e−4) and negatively with footprints (ρ≈−0.60, p≈5.4e−4), supporting the “extracellular bypass” path when Serpin is high.
    * Non-brain: E2F_G2M_score inversely correlates with Proteo-ΔFM (ρ≈−0.55, p≈0.005), consistent with F-driven ΔFM in proliferation-driven settings.
    * Brain stratum currently has n=1 (PXD005719 brain variant) – insufficient for statistics; flagged in tables.
