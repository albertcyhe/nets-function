#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

datasets <- c('GSE184869','GSE125989','GSE43837','GSE14017','GSE14018')

read_scores <- function(ds){
  path <- file.path('data/processed', ds, paste0(ds, '.scores.tsv'))
  fread(path, sep='\t', header=TRUE, data.table=FALSE, check.names=FALSE)
}

sp_cor <- function(x, y){
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 3) return(list(r=NA_real_, p=NA_real_, n=length(x)))
  ct <- suppressWarnings(cor.test(x, y, method='spearman'))
  list(r=unname(ct$estimate), p=ct$p.value, n=length(x))
}

main <- function(){
  out <- data.frame()
  for (ds in datasets){
    sc <- read_scores(ds)
    # singscore F vs M
    r1 <- sp_cor(sc$singscore_F, sc$singscore_M)
    out <- rbind(out, data.frame(dataset=ds, method='singscore', rho=r1$r, p=r1$p, n=r1$n))
    # ssGSEA F vs M
    r2 <- sp_cor(sc$ssgsea_F, sc$ssgsea_M)
    out <- rbind(out, data.frame(dataset=ds, method='ssgsea', rho=r2$r, p=r2$p, n=r2$n))
  }
  out$fdr <- p.adjust(out$p, method='BH')
  dir.create('results/tables', recursive=TRUE, showWarnings=FALSE)
  write.table(out, file='results/tables/fm_correlation_overview.tsv', sep='\t', quote=FALSE, row.names=FALSE)
  message('Wrote: results/tables/fm_correlation_overview.tsv')
}

if (identical(environment(), globalenv())){
  tryCatch(main(), error=function(e){ message('Error: ', e$message); quit(status=1) })
}

