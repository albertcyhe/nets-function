#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GSVA)
})

ensure_pkg <- function(pkg){
  if (!requireNamespace(pkg, quietly=TRUE)) return(FALSE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  TRUE
}

read_expr <- function(path){
  message('Reading expr: ', path)
  df <- read.delim(path, check.names = FALSE)
  rn <- df[[1]]; df <- df[,-1,drop=FALSE]; rownames(df) <- rn
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(as.character(x))))
  as.matrix(df)
}

map_symbols_if_needed <- function(expr){
  rn <- rownames(expr)
  if (length(rn)==0) return(expr)
  if (grepl('^ENSG', rn[1])){
    ok <- FALSE
    try({
      if (!ensure_pkg('AnnotationDbi') || !ensure_pkg('org.Hs.eg.db')) stop('anno missing')
      ens <- vapply(strsplit(rn, '.', fixed = TRUE), function(x) x[[1]], character(1))
      map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=unique(ens), columns=c('SYMBOL'), keytype='ENSEMBL')
      map <- map[!is.na(map$SYMBOL) & map$SYMBOL!='', c('ENSEMBL','SYMBOL')]
      idx <- match(ens, map$ENSEMBL)
      sy <- map$SYMBOL[idx]
      sy[is.na(sy)] <- rn[is.na(sy)]
      rownames(expr) <- toupper(sy)
      expr <- t(apply(expr, 2, function(col){ tapply(col, INDEX = rownames(expr), FUN = function(v) max(v, na.rm=TRUE)) }))
      expr <- t(expr)
      ok <- TRUE
    }, silent=TRUE)
    if (!ok) return(expr)
  }
  rownames(expr) <- toupper(rownames(expr))
  expr
}

collapse_duplicated_rows <- function(expr){
  if (any(duplicated(rownames(expr)))){
    expr <- t(apply(expr, 2, function(col){
      tapply(col, INDEX = rownames(expr), FUN = function(v) max(v, na.rm=TRUE))
    }))
    expr <- t(expr)
  }
  expr
}

get_gene_sets <- function(){
  # Prefer cached msigdbr datasets if available
  if (file.exists('resources/pathways/msig_h.rds')){
    msig_h <- readRDS('resources/pathways/msig_h.rds')
  } else {
    if (!ensure_pkg('msigdbr')) stop('msigdbr not available')
    msig_h <- msigdbr::msigdbr(species = 'Homo sapiens', collection = 'H')
  }
  needed <- c('HALLMARK_TNFA_SIGNALING_VIA_NFKB',
              'HALLMARK_IL6_JAK_STAT3_SIGNALING',
              'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
              'HALLMARK_G2M_CHECKPOINT',
              'HALLMARK_E2F_TARGETS',
              'HALLMARK_ANGIOGENESIS',
              'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')
  h_sets <- split(toupper(msig_h$gene_symbol), msig_h$gs_name)
  h_sets <- h_sets[names(h_sets) %in% needed]

  # p38 MAPK and IL-1 response from Reactome/GO
  msig_c2r <- if (file.exists('resources/pathways/msig_c2_reactome.rds')) readRDS('resources/pathways/msig_c2_reactome.rds') else msigdbr::msigdbr(species='Homo sapiens', collection='C2', subcollection='REACTOME')
  re_list <- split(toupper(msig_c2r$gene_symbol), msig_c2r$gs_name)
  p38_keys <- grep('P38|P-38', names(re_list), ignore.case=TRUE, value=TRUE)
  p38_sets <- re_list[p38_keys]

  msig_c5 <- if (file.exists('resources/pathways/msig_c5_gobp.rds')) readRDS('resources/pathways/msig_c5_gobp.rds') else msigdbr::msigdbr(species='Homo sapiens', collection='C5', subcollection='GO:BP')
  go_list <- split(toupper(msig_c5$gene_symbol), msig_c5$gs_name)
  il1_keys <- grep('INTERLEUKIN_1|IL1', names(go_list), ignore.case=TRUE, value=TRUE)
  il1_sets <- go_list[il1_keys]

  sets <- c(h_sets, p38_sets, il1_sets)
  # de-duplicate genes and drop tiny sets
  sets <- lapply(sets, function(g){ unique(g) })
  sets <- sets[sapply(sets, length) >= 10]
  sets
}

score_pathways <- function(expr, sets){
  # intersect genes
  sets <- lapply(sets, function(g) intersect(g, rownames(expr)))
  sets <- sets[sapply(sets, length) >= 10]
  if (length(sets)==0) return(NULL)
  param <- GSVA::ssgseaParam(exprData = expr, geneSets = sets, normalize = TRUE)
  mat <- GSVA::gsva(param)
  mat
}

main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) stop('Usage: Rscript 05_pathway_scores.R --dataset <ID> --expr <expr.tsv.gz> --outdir <dir>')
  ds <- NULL; expr_path <- NULL; outdir <- NULL
  i <- 1
  while (i <= length(args)){
    k <- args[[i]]
    if (k %in% c('--dataset','--expr','--outdir')){
      v <- args[[i+1]]; i <- i + 1
      if (k=='--dataset') ds <- v
      if (k=='--expr') expr_path <- v
      if (k=='--outdir') outdir <- v
    }
    i <- i + 1
  }
  if (is.null(ds) || is.null(expr_path) || is.null(outdir)) stop('missing args')
  dir.create(outdir, recursive=TRUE, showWarnings = FALSE)

  expr <- read_expr(expr_path)
  expr <- map_symbols_if_needed(expr)
  expr <- collapse_duplicated_rows(expr)
  sets <- get_gene_sets()
  mat <- score_pathways(expr, sets)
  if (is.null(mat)) stop('No pathway scores computed (empty gene set overlap).')
  out <- file.path(outdir, paste0(ds, '.pathways.tsv'))
  tab <- data.frame(pathway = rownames(mat), mat, check.names = FALSE)
  write.table(tab, out, sep='\t', quote=FALSE, row.names=FALSE)
  message('Wrote: ', out)
}

if (identical(environment(), globalenv())){
  tryCatch(main(), error=function(e){ message('Error: ', e$message); quit(status=1) })
}
