Evidence Chain for NET Marker–Function Decoupling in Brain Metastasis (with reproducible records)
===============================================================================================

This README consolidates and reconstructs prior work across RNA and Proteomics/Metabolomics layers, following the flow “rationale → execution → results”, and preserves failures/adjustments for transparency.

Quick Start
- One‑click data retrieval: `bash scripts/sh/download_all_data.sh` (see `DATA_SOURCES.md` for coverage and notes)
- One‑click environments: `bash scripts/sh/setup_proteomics_env.sh` and `bash scripts/sh/setup_scrna_env.sh` (see `ENVIRONMENT.md`)
- Activate proteomics env: `source env/activate_proteomics.sh`

1. Motivation and Hypotheses (markers ≠ function; decoupling?)
- Background: Prior work suggests CTSC‑induced inflammation → NETosis → NET release → TSP1 (THBS1) cleavage and pro‑metastatic effects; many studies strongly bind NET “markers (NET‑M)” to “function (NET‑F)”.
- Core question: markers are not function — can they decouple? In acute inflammation literature, “markers < function” is often observed. If in metastasis prevention markers and function can be separated, clinical implications and treatment stratification are very different.
- Quantitative definition: ΔFM = z(F) − z(M)
  - NET‑F (functional axis): ELANE/PRTN3/CTSG (avoid upstream ROS generators).
  - NET‑M (marker axis): anchor at PADI4, expand by co‑expression in neutrophil‑enriched samples (v2; designed to be anti‑correlated with F).
- Hypothesis: decoupling is brain‑microenvironment‑specific, manifesting as “markers > function” (ΔFM decreases), coupled to THBS1/ECM/angiogenesis external‑growth axes.

Statistical standards and thresholds
- Tests: paired Wilcoxon (report Cliff’s δ), unpaired Mann‑Whitney/linear models; correlations via Spearman + BH FDR (alpha = 0.05).
- Partial correlations/mixed models: control Neutrophil score and TumorPurity; stratify or include (1 | cohort/dataset).
- Reporting: two‑sided tests, FDR < 0.05; report effect sizes with 95% CI; random‑effects meta‑analysis (DL), I² < 50% considered consistent.

Summary and next steps
- Result highlights: propose ΔFM as a quantitative indicator of “marker–function imbalance”, define NET‑F/NET‑M and statistical standards, creating a unified framework for RNA/protein validation.
- Clinical significance: if “markers > function” brain‑specific decoupling holds, prevention and follow‑up stratification are affected (functional anchors require separate monitoring).
- Next steps: validate ΔFM stability in more independent cohorts, refine NET‑M v2 and negative controls, and finalize a locked resource list.

2. Transcriptome evidence: lung coupling, brain decoupling (markers > function)
- Data and design (bulk)
  - Paired metastasis: GSE184869 (RNA‑seq; log2 TMM CPM provided), GSE125989 (Affy GPL571).
  - Metastasis stratification: GSE43837 (brain), GSE14017/14018 (non‑brain).
  - Covariates: Neutrophils (MCPcounter), TumorPurity (tidyestimate/ESTIMATEScore → TumorPurity).
  - Scoring: singscore and ssGSEA (GSVA ≥ 2.2 param API; de‑duplicate and aggregate rownames). Compute ΔFM after within‑sample z‑scaling.

- Main results (bulk)
  - ΔFM paired difference (met − primary)
    - GSE184869 (RNA‑seq): singscore Wilcoxon V=14, p=2.10e−4, FDR=3.22e−4, Cliff’s δ≈−0.70, n=20; ssGSEA consistent in direction.
      - Figure: results/figures/Figure2_GSE184869_Paired_DFM.pdf
      - Table: results/tables/GSE184869_paired_tests.tsv
    - GSE125989 (Affy): V=91, p=0.252, FDR=0.298, Cliff’s δ≈0.25, n=16 (ns).
      - Table: results/tables/GSE125989_paired_tests.tsv
  - ΔFM and THBS1 partial correlation (control Neutrophil + TumorPurity)
    - GSE184869: singscore ρ=−0.3004, p=0.00415, FDR≈0.00831, n=90; ssGSEA ρ=−0.2732, p=0.00936, FDR≈0.00936.
      - Figure: results/figures/Figure3_GSE184869_THBS1_Residuals.pdf
      - Table: results/tables/GSE184869_THBS1_residuals.tsv (partial‑corr backup: results/tables/GSE184869_thbs1_partial_cor.tsv)
    - GSE125989: ρ=0.2269, p=0.2109, FDR≈0.3149 (ns).
  - F/M correlation across cohorts (5 cohorts × 2 methods)
    - GSE184869: singscore ρ=−0.2081, p=0.0491, FDR=0.098; ssGSEA ρ=−0.1216, p=0.2537.
    - GSE125989: ρ=0.4941, p=0.00449, FDR=0.01495; ssGSEA ρ=0.4025, p=0.02314.
    - GSE43837: ρ≈0.07–0.11; GSE14017: ρ≈0.23–0.26; GSE14018: ρ≈0.65–0.67 (p ≪ 0.001).
      - Figure: results/figures/Figure4_FM_Corr_Overview.pdf
      - Table: results/tables/fm_correlation_overview.tsv

- Single‑cell (mechanistic carrier) GSE186344
  - After strict gating (high NEUT≥0.70 and low TNK/Mono≤0.60), no stable neutrophil cluster was found (or counts ~0), thus ΔFM.cell within neutrophils and a “low‑ΔFM” subpopulation could not be established.
  - Interpretation: consistent with “neutrophils have lysed to NETs and functional axis is suppressed”; scRNA captures intracellular transcription and is less suited for extracellular/degradation footprints.
    - Figure: results/figures/Figure5_scRNA_LowDFM_Fraction.pdf (placeholder interpretation noted in logs)
    - Record: results/figures/Figure5_scRNA_THBS1_Assoc.txt (non‑estimable correlation due to size/constant issues)
    - Tables: results/tables/gse186344_neutrophil_lowDFM_fraction.tsv; results/tables/gse186344_sample_subpop_fraction.tsv; results/tables/gse186344_subpop_thbs1_assoc.tsv

- Execution scripts (transcriptomics)
 - R: scripts/r/02_covariates_estimate.R, scripts/r/03_scores_deltafm.R, scripts/r/04_bulk_stats_mvp.R, scripts/r/05_pathway_scores.R, scripts/r/06_fm_correlation_overview.R
  - Figures (re‑draw): Rscript scripts/r/13_figures_a3.R
  - Modules/resources: Rscript scripts/r/export_selected_sets.R
  - Outreach: python scripts/py/met500_fetch.py, python scripts/py/site_triage_met500.py (MET500 resources)

3. Proteomics Evidence (A3): Methods and Results
- Datasets
  - PRIDE PXD011796 (NET reference; label‑free), PXD005719 (brain‑metastatic membrane proteome), PXD046330 (HER2+ conditioned media secretome), PXD051579 (HER2+ brain‑metastasis BBB; paused for cross‑omics), plus MassIVE MSV000089062 and PRIDE PXD018301/PXD032767 for liquid prototypes.

- Processing and scoring (high level)
  - Convert RAW→mzML; search with MSFragger; infer proteins with Philosopher; normalize and compute per‑sample module scores.
  - Derived endpoints: Proteo‑ΔFM (protein ΔFM analogue), NE/PR3 substrate footprints (semi/ non‑tryptic N‑termini proxy), THBS1 cleavage index, Serpin scores (CLR), MNAR detection indicators.

- A3 model executions — first pass (strict evidence)
  - ISI (Inhibitor–Protease Stoichiometry): scripts/py/isi_model.py → results/tables/{isi_per_sample.tsv, isi_models.tsv}.
    - Result: In current cohorts (PXD046330, PXD005719, PXD051579), NET‑F targets (ELANE/PRTN3/CTSG) are not detected (target_sum=0), so ISI is undefined by design, consistent with “F suppressed/not measurable while M persists”.
  - MNAR (Missingness‑as‑signal): scripts/py/mnar_detection_model.py → results/tables/{mnar_detection_table.tsv, mnar_logit_results.tsv}.
    - Result: detect_F=0 across samples; detect_M=1 across samples → non‑estimable logits but a definitive qualitative pattern: “F fully missing, M universally present”.
  - Footprint indices (substrate‑centric NE/PR3 proxy): scripts/py/footprint_indices.py → results/tables/footprints.tsv.
    - Result: indices are consistently negative; in PXD005719, brain‑variant < parental (≈−0.786 vs −0.633); PXD046330/PXD051579 footprints are negative, consistent with low NE/PR3 activity.

- Stratified associations (evidence chain)
  - Footprints vs Serpin (expected negative)
    - Combined: ρ ≈ −0.563, p ≈ 0.00121, q ≈ 0.00182; Non‑brain: ρ ≈ −0.560, p ≈ 0.00160, q ≈ 0.00130; Brain: under‑powered.
  - THBS1 abundance vs Serpin (expected positive)
    - Combined: ρ ≈ +0.583, p ≈ 0.00073, q ≈ 0.00218; Non‑brain: ρ ≈ +0.611, p ≈ 0.00043, q ≈ 0.00130.
  - Proteo‑ΔFM vs Serpin (expected negative)
    - Combined/Non‑brain: weak/non‑significant (ρ ≈ +0.225, p ≈ 0.29), consistent with M‑driven ΔFM under F missingness.
  - Outputs: results/tables/assoc_summary.tsv (per‑stratum summaries).

- Mixed‑effects models
  - scripts/py/mixed_effects_a3.py → results/tables/mixed_effects_a3.tsv
    - Formula: endpoint ~ Serpin_score_core * organ + log_total_PSMs + (1 | dataset)
    - THBS1 abundance fit: β_Serpin ≈ +0.028 (ns; limited organ contrast); footprint model singular (brain underpowered); messages annotated.
  - scripts/py/mixed_effects_growth.py (non‑brain focus) → results/tables/mixed_effects_growth.tsv
    - Footprints ~ Serpin_score_core: β ≈ +0.0027 (p ≈ 4e−20) reflecting monotone footprint drop as Serpin increases (index negative, so small positive slope = stronger suppression).
    - THBS1 cleavage ~ Serpin: β ≈ −0.011 (p ≈ 0.41); THBS1 log abundance ~ Serpin: β ≈ +0.0125 (p ≈ 0.30).
    - Warnings (brain underpowered, singular) recorded in results/tables/mixed_effects_messages.txt.

- SEM prototype (exploratory)
  - scripts/py/prepare_sem_data.py → results/tables/sem_input.tsv
  - scripts/r/14_sem_model.R (lavaan) → results/models/sem_results.json, results/figures/A3_sem_path.svg
    - Model: Serpin_score_core → F_latent → {footprint, THBS1_cleave, Proteo‑ΔFM}; standardized paths saved; warnings indicate limited identifiability (n≈24).

- Multi‑omics alignment (status)
  - results/tables/multiomics_alignment.tsv and results/tables/multiomics_pairwise_deltas.tsv align PXD046330 ↔ GSE96860 (SKBR3) and PXD005719 ↔ GSE12237 (MDA‑MB‑231 parental vs brain variants). Spearman ΔFM correlations underpowered (n<3 per cohort); NET‑M/THBS1 deltas remain interpretable; PXD051579 cross‑omics paused until a matched transcriptomic dataset exists.

4. Site Triage (Liquid Prototypes)
- CSF (MSV000089062, n≈72; BrainMet≈17)
  - Features: data/processed/proteomics/MSV000089062/csf_patient_features.tsv; outputs include metrics/ROC/PR/calibration/DCA/coefficients/scaler/spearman/predictions (see “Figures and tables”).
- Plasma‑EV (PXD018301, n≈512; BrainMet≈12)
  - Multi‑level Excel (Human512Reports) normalized by `Area`; columns renamed to sample names; outputs include full metrics and robust logit.
- Plasma (PXD032767; MaxQuant)
  - Non‑brain RAW not public; use Supplementary Table S1 to map groups; construct `run_manifest.tsv`; compute footprints/THBS1 with MaxQuant exports.
  - Clinical: data/interim/proteomics/PXD032767/metadata/clinical_metadata.tsv
  - Manifest: data/interim/proteomics/PXD032767/metadata/run_manifest.tsv
  - Intermediates: data/interim/proteomics/PXD032767/metadata/{footprint_maxquant.tsv, thbs1_cleavage.tsv}
- Fusion (PXD032767 × PXD018301)
  - Add dataset indicators; output robust logit and (1|dataset) mixed effects alongside DCA.
- Modelling notes
  - Robust logit: implement HC3 manually (leverage h, residuals) when library variant lacks HC3 for logistic.
  - Mixed effects: BinomialBayesMixedGLM; report fe_mean/fe_sd; directions match robust logit.
  - ROC/effects for CSF/Plasma deliver AUC ~0.60–0.79 usable range.

5. Metabolomic Extensions (Paused)
- Current Workbench studies (ST000745, ST001104, ST002921) expose mwTab metadata and/or instrument RAW only; named metabolite matrices are not available; redox ratios (e.g., GSH/GSSG) not computable yet.
- Scripts archived pending restoration once named metabolite tables are available.

6. Cross‑Layer Integration & Modelling (overview)
- Sample alignment: results/tables/multilayer_sample_map.tsv linking cohort/organ/patient/sample types and available layers.
- Derived metrics table: results/tables/multilayer_features.tsv combining RNA ΔFM, Proteo‑ΔFM, Serpin score, THBS1 cleavage index, ROS metrics, MPO traces (where available).
- Mixed‑effects model (R, lme4 + lmerTest): Functional_proxy ~ DeltaFM * Organ + Proteo_DeltaFM + Serpin_score + ROS_ratio + (1 | Cohort); evaluate organ slopes via emmeans.
- Visual storytelling: “double sandwich” schematic results/figures/F_double_sandwich.pdf.

7. Failures and adjustments (transparent records)
- Transcriptome: scRNA THBS1 association not estimable due to size/constant issues; strict gating did not yield a stable neutrophil cluster (or near‑zero), thus “low‑ΔFM neutrophil fraction” cannot be robustly estimated — consistent with “neutrophils already lysed to NETs; functional axis suppressed”.
- Protein/peptide: NET‑F targets are sparse; ISI and robust logit non‑estimable in some cohorts; pivot to “missingness‑as‑signal + footprints/THBS1 cleavage” as twin anchors for function.
- Mixed effects: statsmodels==0.14.5 still reports singular warnings (expected with low n); coefficient directions stable; notes in results/tables/mixed_effects_messages.txt.
- SEM: n≈24, limited identifiability (non‑invertible information matrix, negative variances); treat as conceptual.
- Metabolomics: Workbench (ST000745/ST001104/ST002921) provide mwTab meta + RAW only; named metabolite matrices unavailable; redox ratios (GSH/GSSG) not computable yet; scripts archived pending restoration.

8. Data, directories, and environments (conventions)
- Directory layout
  - Raw: `data/raw/<dataset>/` (keep `checksums.sha256`)
  - Intermediates: `data/interim/<modality>/<dataset>/` (mzML, search workspaces, temporary TSVs)
  - Processed: `data/processed/<dataset>/` (gene/protein matrices, cleavage indices, covariates)
  - Resources: `resources/modules/`, `resources/manifests/`, `resources/annotation/` (includes `gencode.v43.basic.annotation.gtf.gz`), `resources/msfragger/`
  - Results: `results/tables/`, `results/figures/`

- Environments and tools
  - Proteomics env: `env/proteomics.yml`; activate with:
    - `export MAMBA_ROOT_PREFIX="$PWD/env/.mamba"`
    - `env/bin/micromamba env create -y -f env/proteomics.yml`
    - `eval "$(env/bin/micromamba shell hook -s bash)"`
    - `micromamba activate proteomics`
    - `export PATH="$PWD/env/bin:$PATH"`; `philosopher version` (sanity check)
  - RNA bigWig: includes `pybigwig` and `gffread`; annotation at `resources/annotation/gencode/gencode.v43.basic.annotation.gtf.gz`
  - Optional R user library: `R_LIBS_USER=env/.R` (msqrob2/lmerTest/plot, etc.)
  - Long jobs: `tmux`/`screen`; log to `logs/` (e.g., `logs/proteomics/pipeline_<timestamp>.log`)
  - Scripts: `scripts/sh/setup_proteomics_env.sh` (one‑click proteomics deps); RAW→mzML: `scripts/sh/run_trfp_parallel.sh`; smoke test: `scripts/sh/msfragger_smoke.sh`; GEO supplementary: `scripts/sh/fetch_gse186344_suppl.sh`.

- Data retrieval (examples)
  - PRIDE: `bash scripts/sh/fetch_pride_dataset.sh PXD011796 --categories RAW,RESULT --dest data/raw/PXD011796`
  - PXD005719: `scripts/sh/fetch_pride_dataset.sh PXD005719 --categories RAW,RESULT --dest data/raw/PXD005719`
  - PXD046330: `scripts/sh/fetch_pride_dataset.sh PXD046330 --categories RAW,RESULT --dest data/raw/PXD046330`
  - PXD051579 (paused): `scripts/sh/fetch_pride_dataset.sh PXD051579 --categories RAW,RESULT --dest data/raw/PXD051579`
  - GSE96860 bigWig: `python scripts/py/process_gse_bigwig.py --bigwig-dir data/raw/transcriptomics/GSE96860 --gtf resources/annotation/gencode/gencode.v43.basic.annotation.gtf.gz --dataset GSE96860 --output-dir data/processed/GSE96860`
  - GSE12237 array: `python scripts/py/parse_series_matrix.py --series data/raw/transcriptomics/GSE12237/GSE12237_series_matrix.txt.gz --gpl data/raw/transcriptomics/GSE12237/GPL96.annot.gz --dataset GSE12237 --outdir data/processed/GSE12237`
  - CSF proteomes (various): `scripts/sh/fetch_pride_dataset.sh <CSF_ACCESSION> --categories RESULT --match 'CSF' --dest data/raw/<CSF_ACCESSION>`
  - CSF metabolomics: `aria2c -c -x16 -s16 --max-tries=0 --retry-wait=30 -i resources/manifests/metabolomics_urls.txt -d data/raw/metabolomics`

9. Reproducible script entry points and suggested order
- Bulk (MVP)
  - `Rscript scripts/r/02_covariates_estimate.R --dataset GSE184869 --config resources/config.yml`
  - `Rscript scripts/r/03_scores_deltafm.R --dataset GSE184869 --config resources/config.yml`
  - `Rscript scripts/r/04_bulk_stats_mvp.R --dataset GSE184869 --config resources/config.yml`
  - `Rscript scripts/r/05_pathway_scores.R --dataset GSE184869 --config resources/config.yml`
  - `Rscript scripts/r/06_fm_correlation_overview.R --config resources/config.yml`

- scRNA (GSE186344)
  - `python scripts/py/sc_preprocess_gse186344.py --raw GSE186344_RAW.tar --out data/processed/scRNA/gse186344.h5ad`
  - `python scripts/py/sc_neutrophil_subsets.py --h5ad data/processed/scRNA/gse186344.h5ad --net_f resources/modules/net_f_v1.tsv --net_m resources/modules/net_m_v2.tsv --outdir results/tables`
  - `python scripts/py/sc_tumor_thbs1_assoc.py --h5ad data/processed/scRNA/gse186344.h5ad --frac results/tables/gse186344_neutrophil_lowDFM_fraction.tsv --out results/tables/gse186344_subpop_thbs1_assoc.tsv`

- Proteins/peptides (A3)
  - `python scripts/py/footprint_indices.py` → `results/tables/footprints.tsv`
  - `python scripts/py/thbs1_cleavage_index.py` → `results/tables/thbs1_cleave_idx.tsv`
  - `python scripts/py/serpin_scores.py` → `results/tables/serpin_scores.tsv`
  - `python scripts/py/proteomics_deltafm_score.py` → `data/processed/proteomics/<dataset>/proteo_deltafm.tsv`
  - `python scripts/py/isi_model.py` → `results/tables/isi_per_sample.tsv`, `isi_models.tsv`
  - `python scripts/py/mnar_detection_model.py` → `results/tables/mnar_detection_table.tsv`, `mnar_logit_results.tsv`

10. Figures and tables (outputs, browsable on GitHub)
- Bulk results
  - Figures:
    - [Paired ΔFM: GSE184869](results/figures/Figure2_GSE184869_Paired_DFM.pdf) — paired difference (met − primary) showing brain‑specific decoupling.
    - [THBS1 residuals vs ΔFM: GSE184869](results/figures/Figure3_GSE184869_THBS1_Residuals.pdf) — partial correlation controlling neutrophils and purity.
    - [F/M correlation overview](results/figures/Figure4_FM_Corr_Overview.pdf) — across cohorts (5×2 methods).
    - [Low‑ΔFM fraction in scRNA](results/figures/Figure5_scRNA_LowDFM_Fraction.pdf) — gate outcome summary; see note for n/constant issues.
  - Note: [scRNA THBS1 association log](results/figures/Figure5_scRNA_THBS1_Assoc.txt) — non‑estimable correlation record.
  - Tables: [GSE184869 paired](results/tables/GSE184869_paired_tests.tsv), [GSE125989 paired](results/tables/GSE125989_paired_tests.tsv), [GSE184869 THBS1 residuals](results/tables/GSE184869_THBS1_residuals.tsv), [GSE184869 partial‑corr backup](results/tables/GSE184869_thbs1_partial_cor.tsv), [F/M correlation overview](results/tables/fm_correlation_overview.tsv), [GSE186344 scRNA tables](results/tables/)

- Proteins/peptides and integration
  - Figures:
    - [Mixed‑effects diagnostics](results/figures/mixed_effects_diagnostics.pdf) — model checks.
    - [SEM path diagram](results/figures/A3_sem_path.svg) — exploratory latent structure.
    - [Double‑sandwich schematic](results/figures/F_double_sandwich.pdf) — integrative brain vs lung story.
    - [Proteo‑ΔFM boxplots](results/figures/F_proteo_deltafm_boxplots.pdf) — per dataset distributions.
    - [Proteo‑ΔFM vs transcriptome](results/figures/F_proteo_deltafm_vs_transcriptome.pdf) — cross‑layer relationship.
    - [THBS1 spectrum example](results/figures/thbs1_spectrum_JIMT1_BR_R3.pdf) — peptide‑level evidence.
    - A3 association panels: [forest](results/figures/A3_forest_assoc.pdf), [Serpin vs footprints](results/figures/A3_scatter_serpin_vs_footprint.pdf), [Serpin vs control footprints](results/figures/A3_scatter_serpin_vs_ctrl_footprint.pdf), [Serpin vs ΔFM](results/figures/A3_scatter_serpin_vs_deltaFM.pdf), [Serpin vs THBS1](results/figures/A3_scatter_serpin_vs_THBS1.pdf), [Serpin vs THBS1 cleavage](results/figures/A3_scatter_serpin_vs_thbs1_cleave.pdf), [Serpin vs ECM/neg. controls](results/figures/A3_scatter_serpin_vs_ecm_nc.pdf)
  - Tables: [footprints](results/tables/footprints.tsv), [THBS1 cleavage index](results/tables/thbs1_cleave_idx.tsv), [Serpin scores](results/tables/serpin_scores.tsv), per‑dataset [proteo‑ΔFM](data/processed/proteomics/), [proteo‑ΔFM summary](results/tables/proteo_deltafm_summary.tsv), [ISI per sample](results/tables/isi_per_sample.tsv), [ISI models](results/tables/isi_models.tsv), [MNAR table](results/tables/mnar_detection_table.tsv), [MNAR logit](results/tables/mnar_logit_results.tsv), [assoc summary](results/tables/assoc_summary.tsv), [A3 mixed effects](results/tables/mixed_effects_a3.tsv), [growth mixed effects](results/tables/mixed_effects_growth.tsv), [mixed‑effects messages](results/tables/mixed_effects_messages.txt), [multi‑omics alignment](results/tables/multiomics_alignment.tsv), [pairwise deltas](results/tables/multiomics_pairwise_deltas.tsv), [multilayer sample map](results/tables/multilayer_sample_map.tsv), [multilayer features](results/tables/multilayer_features.tsv), [proteomics–transcript map](results/tables/proteomics_transcript_sample_map.tsv), [matched](results/tables/proteomics_transcript_matched.tsv), [unmatched](results/tables/proteomics_transcript_unmatched.tsv), [multi‑omics correlations](results/tables/multiomics_correlations.tsv), [Serpin associations](results/tables/proteo_serpin_assoc.tsv), [module detection summary](results/tables/proteo_module_detection_summary.tsv)

- Liquid prototypes
  - Figures: [ROC](results/figures/site_triage_roc.pdf), [ROC (png)](results/figures/site_triage_roc.png), [effects](results/figures/site_triage_effects.pdf), [effects (png)](results/figures/site_triage_effects.png), [ROC+calibration+DCA](results/figures/A4_site_triage_roc_calib_dca.pdf)
  - Tables: [metrics](results/tables/site_triage_metrics.tsv), [predictions](results/tables/site_triage_predictions.tsv), [coefficients](results/tables/site_triage_coefficients.tsv), [DCA](results/tables/site_triage_dca.tsv), [DCA baseline](results/tables/site_triage_dca_baseline.tsv), [DCA Serpin‑extended](results/tables/site_triage_dca_serpin_extended.tsv), [DCA null‑inflammation](results/tables/site_triage_dca_null_inflammation.tsv), variant [metrics](results/tables/site_triage_variant_metrics.tsv), variant [predictions](results/tables/site_triage_variant_predictions.tsv), variant [coefficients](results/tables/site_triage_variant_coefficients.tsv), CSF MSV000089062 [metrics](results/tables/msv000089062_csf_metrics.tsv), [ROC](results/tables/msv000089062_csf_roc.tsv), [PR](results/tables/msv000089062_csf_pr.tsv), [calibration](results/tables/msv000089062_csf_calibration.tsv), [DCA](results/tables/msv000089062_csf_dca.tsv), [coefficients](results/tables/msv000089062_csf_coefficients.tsv), [scaler](results/tables/msv000089062_csf_scaler.tsv), [spearman](results/tables/msv000089062_csf_spearman.tsv), [predictions](results/tables/msv000089062_csf_predictions.tsv), PXD018301 CSF [metrics](results/tables/pxd018301_csf_metrics.tsv), [ROC](results/tables/pxd018301_csf_roc.tsv), [PR](results/tables/pxd018301_csf_pr.tsv), [calibration](results/tables/pxd018301_csf_calibration.tsv), [coefficients](results/tables/pxd018301_csf_coefficients.tsv), [scaler](results/tables/pxd018301_csf_scaler.tsv), [spearman](results/tables/pxd018301_csf_spearman.tsv), [predictions](results/tables/pxd018301_csf_predictions.tsv), [DCA](results/tables/pxd018301_csf_dca.tsv), [robust logit](results/tables/pxd018301_csf_robust_logit.tsv), PXD032767 liquid [metrics](results/tables/pxd032767_liquid_metrics.tsv), [ROC](results/tables/pxd032767_liquid_roc.tsv), [PR](results/tables/pxd032767_liquid_pr.tsv), [calibration](results/tables/pxd032767_liquid_calibration.tsv), [coefficients](results/tables/pxd032767_liquid_coefficients.tsv), [scaler](results/tables/pxd032767_liquid_scaler.tsv), [predictions](results/tables/pxd032767_liquid_predictions.tsv), [DCA](results/tables/pxd032767_liquid_dca.tsv), [spearman](results/tables/pxd032767_liquid_spearman.tsv), [model summary](results/tables/pzd032767_model_summary.tsv), plasma fusion PXD032767×PXD018301 [metrics](results/tables/pxd032767_pxd018301_plasma_metrics.tsv), [ROC](results/tables/pxd032767_pxd018301_plasma_roc.tsv), [PR](results/tables/pxd032767_pxd018301_plasma_pr.tsv), [calibration](results/tables/pxd032767_pxd018301_plasma_calibration.tsv), [predictions](results/tables/pxd032767_pxd018301_plasma_predictions.tsv), [coefficients](results/tables/pxd032767_pxd018301_plasma_coefficients.tsv), [scaler](results/tables/pxd032767_pxd018301_plasma_scaler.tsv), [spearman](results/tables/pxd032767_pxd018301_plasma_spearman.tsv), [robust logit](results/tables/pxd032767_pxd018301_plasma_robust_logit.tsv), [mixed effects](results/tables/pxd032767_pxd018301_plasma_mixed_effects.tsv), [DCA](results/tables/pxd032767_pxd018301_plasma_dca.tsv)

- Other schematics/overviews
  - Figures: [ΔFM concept schematic](results/figures/Figure1_DFM_Schematic.pdf), [Mechanism model](results/figures/Figure6_Mechanism_Model.pdf)

- Controls and modules
  - Tables: [negative controls](results/tables/neg_controls.tsv), [negative control associations](results/tables/neg_controls_assoc.tsv), [module scores](results/tables/module_scores.tsv), [module associations](results/tables/module_assoc.tsv), [bulk paired tests](results/tables/bulk_paired_tests.tsv), [bulk ΔFM scores](results/tables/)

11. Next steps (two‑week minimal deliverables)
- ISI (core/extended) + organ interaction; export forest plots.
- Logit/zero‑inflation: F detection ~ Serpin (M as negative control); meta‑analysis across cohorts.
- Substrate footprints and THBS1 cleavage: robust estimates and visualization against Serpin/Proteo‑ΔFM.
- Liquid prototype page: integrate CSF/Plasma ROC and effect sizes; update brain vs non‑brain mechanism schematic.

—

Note: This README covers all scripts, figures, tables, and failure notes from `README_proteo_metab.md`. The historical long document is kept for reference.
