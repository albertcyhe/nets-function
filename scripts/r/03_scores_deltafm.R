#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GSVA)
  library(singscore)
})

ensure_pkg <- function(pkg, bioc=FALSE){
  if (!requireNamespace(pkg, quietly=TRUE)){
    if (bioc){
      if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(pkg, update=FALSE, ask=FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

read_expr <- function(path){
  message("Reading expr: ", path)
  df <- read.delim(path, check.names = FALSE)
  rn <- df[[1]]
  df <- df[,-1, drop=FALSE]
  rownames(df) <- rn
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(as.character(x))))
  as.matrix(df)
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

map_symbols_if_needed <- function(expr){
  rn <- rownames(expr)
  if (length(rn) == 0) return(expr)
  # If looks like ENSEMBL, map to SYMBOL
  if (grepl("^ENSG", rn[1])){
    ok <- FALSE
    try({
      ensure_pkg("AnnotationDbi", bioc=TRUE)
      ensure_pkg("org.Hs.eg.db", bioc=TRUE)
      ens <- gsub("\\..*$", "", rn)
      map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=unique(ens), columns=c("SYMBOL"), keytype="ENSEMBL")
      map <- map[!is.na(map$SYMBOL) & map$SYMBOL!="", c("ENSEMBL","SYMBOL")]
      idx <- match(ens, map$ENSEMBL)
      sy <- map$SYMBOL[idx]
      sy[is.na(sy)] <- rn[is.na(sy)]
      rownames(expr) <- toupper(sy)
      # collapse duplicates by max
      expr <- t(apply(expr, 2, function(col){ tapply(col, INDEX = rownames(expr), FUN = function(v) max(v, na.rm=TRUE)) }))
      expr <- t(expr)
      ok <- TRUE
    }, silent = TRUE)
    if (!ok) return(expr)
  }
  rownames(expr) <- toupper(rownames(expr))
  expr
}

read_genes <- function(path){
  g <- readLines(path)
  g <- toupper(trimws(g))
  g <- g[g != ""]
  unique(g)
}

zs <- function(x){ s <- sd(x, na.rm=TRUE); m <- mean(x, na.rm=TRUE); ifelse(s>0, (x-m)/s, 0) }

score_modules <- function(expr, f_genes, m_genes){
  # singscore
  common_f <- intersect(rownames(expr), f_genes)
  common_m <- intersect(rownames(expr), m_genes)
  ss_f <- rep(NA_real_, ncol(expr)); ss_m <- rep(NA_real_, ncol(expr))
  names(ss_f) <- colnames(expr); names(ss_m) <- colnames(expr)
  rankData <- tryCatch({ singscore::rankGenes(expr) }, error=function(e) NULL)
  if (length(common_f) >= 1){
    sc <- if (!is.null(rankData)) singscore::simpleScore(rankData, upSet = common_f, centerScore=TRUE) else NULL
    ss_f <- sc$TotalScore
  }
  if (length(common_m) >= 1){
    sc <- if (!is.null(rankData)) singscore::simpleScore(rankData, upSet = common_m, centerScore=TRUE) else NULL
    ss_m <- sc$TotalScore
  }
  # ssGSEA via GSVA (new API in GSVA >= 2: param object)
  gs_list <- list(NET_F = common_f, NET_M = common_m)
  gs_list <- gs_list[sapply(gs_list, length) > 0]
  ssgsea_f <- rep(NA_real_, ncol(expr)); ssgsea_m <- rep(NA_real_, ncol(expr))
  names(ssgsea_f) <- colnames(expr); names(ssgsea_m) <- colnames(expr)
  if (length(gs_list) > 0){
    ssg <- tryCatch({
      param <- GSVA::ssgseaParam(exprData = expr, geneSets = gs_list, normalize = TRUE, checkNA = "auto")
      GSVA::gsva(param)
    }, error=function(e) { message("GSVA/ssGSEA failed: ", e$message); NULL })
    if (!is.null(ssg)){
      mat <- tryCatch({
        if (is.matrix(ssg)) ssg else GSVA::gsvaScores(ssg)
      }, error=function(e) NULL)
      if (!is.null(mat)){
        if ("NET_F" %in% rownames(mat)) ssgsea_f <- as.numeric(mat["NET_F", colnames(expr)])
        if ("NET_M" %in% rownames(mat)) ssgsea_m <- as.numeric(mat["NET_M", colnames(expr)])
      }
    }
  }
  data.frame(sample_id = colnames(expr),
             singscore_F = as.numeric(ss_f),
             singscore_M = as.numeric(ss_m),
             ssgsea_F = as.numeric(ssgsea_f),
             ssgsea_M = as.numeric(ssgsea_m),
             stringsAsFactors = FALSE)
}

main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(dataset=NULL, expr=NULL, outdir=NULL, f_genes="resources/modules/net_f_v1.tsv", m_genes="resources/modules/net_m_v1.tsv")
  i <- 1
  while (i <= length(args)){
    k <- args[[i]]
    if (k %in% c("--dataset","--expr","--outdir","--f","--m")){
      v <- args[[i+1]]; i <- i + 1
      if (k == "--dataset") opt$dataset <- v
      if (k == "--expr") opt$expr <- v
      if (k == "--outdir") opt$outdir <- v
      if (k == "--f") opt$f_genes <- v
      if (k == "--m") opt$m_genes <- v
    }
    i <- i + 1
  }
  if (is.null(opt$dataset) || is.null(opt$expr) || is.null(opt$outdir)){
    stop("Usage: Rscript 03_scores_deltafm.R --dataset <ID> --expr <expr.tsv.gz> --outdir <dir> [--f <net_f.tsv>] [--m <net_m.tsv>]")
  }
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

  expr <- read_expr(opt$expr)
  expr <- map_symbols_if_needed(expr)
  f_genes <- read_genes(opt$f_genes)
  m_genes <- read_genes(opt$m_genes)
  expr <- collapse_duplicated_rows(expr)

  scores <- score_modules(expr, f_genes, m_genes)
  # Z-standardize and compute ΔFM per method
  scores$deltaFM_singscore <- zs(scores$singscore_F) - zs(scores$singscore_M)
  scores$deltaFM_ssgsea   <- zs(scores$ssgsea_F)   - zs(scores$ssgsea_M)

  out_scores <- file.path(opt$outdir, paste0(opt$dataset, ".scores.tsv"))
  write.table(scores, out_scores, sep="\t", quote=FALSE, row.names=FALSE)
  message("Wrote: ", out_scores)
}

if (identical(environment(), globalenv())){
  tryCatch(main(), error=function(e){
    message("Error: ", e$message)
    quit(status=1)
  })
}
