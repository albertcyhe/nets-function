Environments and One-Click Setup
================================

This repository uses two conda/micromamba environments:
- `proteomics` — Python + R toolchain for proteomics, pathway scoring, statistics, and utilities.
- `nets-scrna` — Python stack for single‑cell analyses.

One‑Click Setup
- Proteomics toolchain: `bash scripts/sh/setup_proteomics_env.sh`
  - Installs micromamba locally under `env/` (if missing)
  - Creates/updates the `proteomics` env from `env/proteomics.yml`
  - Downloads and wires `philosopher` into `env/bin/`
  - Prints the activation command
- scRNA toolchain: `bash scripts/sh/setup_scrna_env.sh`
  - Reuses/installs micromamba under `env/`
  - Creates/updates the `nets-scrna` env from `env/environment-scrna.yml`

Manual Activation
- Quick activate (proteomics):
  - `source env/activate_proteomics.sh`
- Or via micromamba directly:
  - `source <(env/bin/micromamba shell hook --shell bash)`
  - `micromamba activate proteomics` (or `micromamba activate nets-scrna`)

Verification
- `philosopher version` (after `proteomics` activation)
- `python -c "import numpy, pandas; print(numpy.__version__, pandas.__version__)"`
- `R -q -e 'sessionInfo()'` (optional; validates R stack in `proteomics`)

Primary Specs
- `env/proteomics.yml` includes: Python 3.10, numpy/pandas/scipy/statsmodels, pyteomics/pymzml, ProteoWizard, MSFragger, OpenMS, OpenJDK 17, R 4.3 + ggplot2/dplyr/readr/tidyr/broom, Bioconductor `msqrob2`, rpy2, jq/yq, pybigwig, gffread.
- `env/environment-scrna.yml` includes: scanpy/anndata stack, harmonypy, scrublet, DoubletDetection, and plotting/scikit-learn dependencies.

Tips
- For long runs, use `tmux`/`screen` and capture logs in `logs/`.
- To freeze and share exact versions, export `micromamba env export -n <env> > env/export_<env>.yml`.

