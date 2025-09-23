#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(lavaan)
})

args <- list(
  input = 'results/tables/sem_input.tsv',
  out_json = 'results/models/sem_results.json',
  out_svg = 'results/figures/A3_sem_path.svg'
)

dat <- read_tsv(args$input, show_col_types = FALSE)

if(nrow(dat) < 20){
  warning('SEM data limited; results exploratory only.')
}

# z-standardise predictors for interpretability
scale_cols <- c('Serpin_score_core','footprint_index','THBS1_cleave_idx','Proteo_DeltaFM')
for (col in scale_cols){
  if(col %in% names(dat)){
    dat[[paste0(col,'_z')]] <- scale(dat[[col]])[,1]
  }
}

model <- '
  F_latent =~ footprint_index_z + THBS1_cleave_idx_z + Proteo_DeltaFM_z
  Proteo_DeltaFM_z ~ b1*F_latent
  F_latent ~ a1*Serpin_score_core_z
  # indirect effect
  ind_effect := a1 * b1
'

fit <- sem(model, data = dat, estimator='MLR', missing='fiml')

summary_fit <- parameterEstimates(fit, standardized=TRUE)
fit_measures <- fitMeasures(fit, c('cfi','tli','rmsea','srmr'))

dir.create(dirname(args$out_json), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(args$out_svg), showWarnings = FALSE, recursive = TRUE)

library(jsonlite)
write_json(list(parameters=summary_fit, fit=fit_measures), path=args$out_json, pretty=TRUE, auto_unbox=TRUE)

library(semPlot)
semPaths(fit, what='std', layout='tree', edge.label.cex=1.0, sizeMan=6, sizeLat=7,
         mar=c(4,4,4,4), esize=1.2, residuals=FALSE, intercepts=FALSE,
         label.prop=0.7, pastel=TRUE)
dev.copy(svg, filename=args$out_svg, width=6, height=4)
dev.off()

message('SEM results written to ', args$out_json, ' and ', args$out_svg)

