#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

read_expr <- function(path){
  df <- fread(path, sep='\t', header=TRUE, data.table=FALSE, check.names=FALSE)
  rn <- df[[1]]; df <- df[,-1,drop=FALSE]; rownames(df) <- rn
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(as.character(x))))
  as.matrix(df)
}

read_vec <- function(path){
  x <- readLines(path); x <- toupper(trimws(x)); unique(x[x!=''])
}

zs <- function(x){ s <- sd(x, na.rm=TRUE); m <- mean(x, na.rm=TRUE); ifelse(s>0,(x-m)/s,0) }

main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) stop('Usage: Rscript 01_build_net_m_v2.R <dataset> <expr.tsv.gz> <covariates.tsv> [out_genes resources/modules/net_m_v2.tsv]')
  ds <- args[[1]]; expr_path <- args[[2]]; cov_path <- args[[3]]
  out_genes <- ifelse(length(args)>=4, args[[4]], 'resources/modules/net_m_v2.tsv')
  f_genes <- read_vec('resources/modules/net_f_v1.tsv')
  expr <- read_expr(expr_path)
  cov <- fread(cov_path, sep='\t', header=TRUE, data.table=FALSE, check.names=FALSE)
  rownames(cov) <- cov$sample_id
  common <- intersect(colnames(expr), cov$sample_id)
  expr <- expr[, common, drop=FALSE]
  cov <- cov[common,,drop=FALSE]
  rownames(expr) <- toupper(rownames(expr))

  if (!('PADI4' %in% rownames(expr))) stop('PADI4 not found in expression; choose a dataset where it exists (e.g., GSE125989).')

  # neutrophil-enriched subset (top quartile by Neutrophil_proxy)
  if (!('Neutrophil_proxy' %in% colnames(cov))) stop('Neutrophil_proxy not found in covariates.')
  thr <- quantile(cov$Neutrophil_proxy, 0.75, na.rm=TRUE)
  keep <- rownames(cov)[cov$Neutrophil_proxy >= thr & is.finite(cov$Neutrophil_proxy)]
  if (length(keep) < 8) keep <- rownames(cov)[order(cov$Neutrophil_proxy, decreasing=TRUE)][seq_len(min(10, nrow(cov)))]
  X <- expr[, keep, drop=FALSE]

  target <- as.numeric(X['PADI4',])
  # compute Spearman correlation with all genes
  rho <- apply(X, 1, function(v){ suppressWarnings(suppressWarnings(suppressMessages(cor(v, target, method='spearman', use='pairwise.complete.obs')))) })
  rho[is.na(rho)] <- 0
  # filter: positive correlation, exclude NET-F genes
  rho <- rho[setdiff(names(rho), f_genes)]
  cand <- sort(rho, decreasing=TRUE)
  # threshold: rho >= 0.3 (relax if too few)
  sel <- names(cand)[cand >= 0.3]
  if (length(sel) < 5) sel <- names(cand)[seq_len(min(10, length(cand)))]
  # sanitize gene symbols (split composite names, remove LOC*)
  clean_sym <- function(x){
    x <- toupper(x)
    x <- sub("[;|/\\s].*$", "", x)
    x
  }
  sel <- clean_sym(sel)
  sel <- sel[sel != '' & !grepl('^LOC', sel)]
  sel <- unique(c('PADI4', sel))
  sel <- setdiff(sel, f_genes)
  # cap to 5-10 genes (prefer 10 if available)
  if (length(sel) > 10) sel <- sel[1:10]
  # write module
  dir.create(dirname(out_genes), recursive=TRUE, showWarnings=FALSE)
  writeLines(sel, out_genes)
  message('Wrote NET-M.v2 genes to: ', out_genes)
  # save correlation table
  out_tab <- file.path(dirname(out_genes), paste0(tolower(ds), '_net_m_v2_cor.tsv'))
  tab <- data.frame(GENE=names(cand), RHO=as.numeric(cand), stringsAsFactors = FALSE)
  write.table(tab, out_tab, sep='\t', quote=FALSE, row.names=FALSE)
  message('Wrote correlations to: ', out_tab)
}

if (identical(environment(), globalenv())){
  tryCatch(main(), error=function(e){ message('Error: ', e$message); quit(status=1) })
}
