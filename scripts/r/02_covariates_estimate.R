#!/usr/bin/env Rscript

ensure_pkg <- function(pkg){
  if (!requireNamespace(pkg, quietly=TRUE)) return(FALSE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  TRUE
}

read_expr <- function(path){
  message("Reading expr: ", path)
  df <- read.delim(path, check.names = FALSE)
  rn <- df[[1]]
  df <- df[,-1, drop=FALSE]
  rownames(df) <- rn
  # coerce to numeric
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(as.character(x))))
  df <- as.matrix(df)
  return(df)
}

map_symbols <- function(expr){
  rn <- rownames(expr)
  if (!all(grepl("^[A-Za-z0-9_.-]+$", rn))) return(expr)
  if (!grepl("^ENSG", rn[1])) return(expr)
  # Try to map Ensembl -> SYMBOL via org.Hs.eg.db
  sym_mat <- expr
  ok <- FALSE
  try({
    if (!ensure_pkg("AnnotationDbi") || !ensure_pkg("org.Hs.eg.db")) stop("anno not available")
    ens <- gsub("\\..*$", "", rn)
    keytypes <- AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)
    if ("ENSEMBL" %in% keytypes){
      map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=unique(ens), columns=c("SYMBOL"), keytype="ENSEMBL")
      map <- map[!is.na(map$SYMBOL) & map$SYMBOL!="", c("ENSEMBL","SYMBOL")]
      map$ENSEMBL <- make.unique(map$ENSEMBL)
      df <- data.frame(ENSEMBL=ens, row=seq_along(ens))
      df$ENSEMBL <- make.unique(df$ENSEMBL)
      df <- merge(df, map, by="ENSEMBL", all.x=TRUE, sort=FALSE)
      sy <- df$SYMBOL
      sy[is.na(sy)] <- rn[is.na(sy)]
      rownames(sym_mat) <- toupper(sy)
      # collapse duplicates by max
      sym_mat <- t(apply(sym_mat, 2, function(col){
        tapply(col, INDEX = rownames(sym_mat), FUN = function(v) max(v, na.rm=TRUE))
      }))
      sym_mat <- t(sym_mat)
      ok <- TRUE
    }
  }, silent = TRUE)
  if (ok) return(sym_mat)
  expr
}

run_mcp <- function(expr){
  res <- NULL
  try({
    if (!ensure_pkg("MCPcounter")) stop("MCPcounter not available")
    # Expect gene symbols in rows
    res <- MCPcounter::MCPcounter.estimate(expr, featuresType = "HUGO_symbols")
  }, silent = TRUE)
  return(res)
}

write_gct <- function(mat, file){
  mat <- cbind(NAME = rownames(mat), Description = rownames(mat), mat)
  con <- file(file, open="wt")
  on.exit(close(con))
  writeLines("#1.2", con)
  writeLines(paste(nrow(mat), ncol(mat)-2, sep="\t"), con)
  write.table(mat, con, sep="\t", quote=FALSE, row.names=FALSE)
}

run_estimate_legacy <- function(expr){
  out <- NULL
  tmpdir <- tempdir()
  in_gct <- file.path(tmpdir, "input.gct")
  filt_gct <- file.path(tmpdir, "input_filtered.gct")
  out_scores <- file.path(tmpdir, "estimate_scores.gct")
  try({
    if (!ensure_pkg("estimate")) stop("estimate not available")
    # estimate expects integers or log-expression; we pass as-is
    write_gct(expr, in_gct)
    estimate::filterCommonGenes(in_gct, filt_gct, id = "GeneSymbol")
    # platform string per docs: "affymetrix" or "agilent" or "illumina"
    estimate::estimateScore(filt_gct, out_scores, platform = "affymetrix")
    # read scores gct
    g <- readLines(out_scores)
    g <- g[-c(1,2)]
    tab <- read.delim(textConnection(paste(g, collapse="\n")), check.names=FALSE)
    rownames(tab) <- tab[[1]]
    tab <- tab[,-c(1,2), drop=FALSE]
    out <- as.data.frame(t(tab))
  }, silent = TRUE)
  return(out)
}

run_tidyestimate <- function(expr, is_affy){
  out <- NULL
  ok <- FALSE
  try({
    if (!ensure_pkg('tidyestimate')) stop('tidyestimate not available')
    # build tibble with first column as hgnc_symbol
    df <- data.frame(hgnc_symbol = rownames(expr), expr, check.names = FALSE)
    # keep only common genes for ESTIMATE
    df2 <- tidyestimate::filter_common_genes(df, id='hgnc_symbol', tidy=TRUE, tell_missing=FALSE)
    sc <- tidyestimate::estimate_score(df2, is_affymetrix = is_affy)
    # rename columns
    colnames(sc) <- sub('^stromal$', 'StromalScore', colnames(sc))
    colnames(sc) <- sub('^immune$', 'ImmuneScore', colnames(sc))
    colnames(sc) <- sub('^estimate$', 'ESTIMATEScore', colnames(sc))
    colnames(sc) <- sub('^purity$', 'TumorPurity', colnames(sc))
    rownames(sc) <- sc$sample
    sc$sample <- NULL
    out <- sc
    ok <- TRUE
  }, silent=TRUE)
  if (!ok) return(NULL)
  out
}

zscore_cols <- function(v){
  m <- mean(v, na.rm=TRUE); s <- sd(v, na.rm=TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(v)))
  (v - m) / s
}

marker_score <- function(expr, genes){
  g <- intersect(rownames(expr), toupper(genes))
  if (length(g) == 0) return(rep(NA_real_, ncol(expr)))
  sub <- expr[g,, drop=FALSE]
  # z-score per gene across samples
  sub <- t(apply(sub, 1, zscore_cols))
  sc <- colMeans(sub, na.rm=TRUE)
  as.numeric(sc)
}

parse_cli <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(dataset = NULL, expr = NULL, outdir = NULL)
  i <- 1
  while (i <= length(args)){
    key <- args[[i]]
    if (key %in% c("--dataset","--expr","--outdir")){
      if (i == length(args)) stop(paste("Missing value for", key))
      val <- args[[i+1]]; i <- i + 1
      if (key == "--dataset") opt$dataset <- val
      if (key == "--expr") opt$expr <- val
      if (key == "--outdir") opt$outdir <- val
    }
    i <- i + 1
  }
  if (is.null(opt$dataset) || is.null(opt$expr) || is.null(opt$outdir)){
    stop("Usage: Rscript 02_covariates_estimate.R --dataset <ID> --expr <expr.tsv.gz> --outdir <dir>")
  }
  opt
}

main <- function(){
  opt <- parse_cli()
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

  expr <- read_expr(opt$expr)
  expr <- map_symbols(expr)

  # MCP-counter
  mcp <- run_mcp(expr)
  if (!is.null(mcp)){
    mcp <- as.data.frame(t(mcp))
  }

  # ESTIMATE/tidyestimate
  # try to infer platform to decide purity computation
  plat_file <- file.path(opt$outdir, paste0(opt$dataset, '.platform.tsv'))
  is_affy <- FALSE
  if (file.exists(plat_file)){
    pf <- tryCatch(read.delim(plat_file, check.names = FALSE), error=function(e) NULL)
    if (!is.null(pf)){
      # heuristic: GEO GPL implies microarray (Affymetrix in our sets); RNA-seq annotated differently
      if ('Series_platform_id' %in% colnames(pf)){
        is_affy <- TRUE
      } else if ('platform' %in% colnames(pf)){
        is_affy <- grepl('affy|hgu|u133', tolower(pf$platform[1]))
      }
    }
  }
  est <- run_tidyestimate(expr, is_affy)
  if (is.null(est)){
    est <- run_estimate_legacy(expr)
  }

  # Merge
  sid <- colnames(expr)
  cov <- data.frame(sample_id = sid, stringsAsFactors = FALSE)
  if (!is.null(mcp)){
    mcp$sample_id <- rownames(mcp)
    cov <- merge(cov, mcp, by = "sample_id", all.x = TRUE)
    # normalize neutrophil column name if present
    nc <- grep('neutro', tolower(colnames(cov)), value=TRUE)
    if (length(nc)>0) cov$Neutrophil <- cov[[nc[1]]]
  }
  if (!is.null(est)){
    est$sample_id <- rownames(est)
    cov <- merge(cov, est, by = "sample_id", all.x = TRUE)
    # If TumorPurity missing but ESTIMATEScore present, derive purity via ESTIMATE formula
    if (!('TumorPurity' %in% colnames(cov)) && ('ESTIMATEScore' %in% colnames(cov))){
      cov$TumorPurity <- cos(0.6049872018 + 0.0001467884 * cov$ESTIMATEScore)
    }
  }

  # Fallback proxies if methods unavailable
  if (is.null(mcp) || !any(grepl("Neutrophil", colnames(cov)))){
    neut_genes <- c("CSF3R","FCGR3B","CXCR2","S100A8","S100A9","ELANE","MPO","LCN2","CEACAM8","MMP8","MMP9")
    ns <- marker_score(expr, neut_genes)
    cov$Neutrophil_proxy <- ns
  }
  if (is.null(est) || !all(c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity") %in% colnames(cov))){
    imm_genes <- c("PTPRC","LCK","CD3D","CD3E","CD4","CD8A","ITGAM","HLA-DRA","HLA-DRB1","CCR7")
    str_genes <- c("COL1A1","COL1A2","COL3A1","DCN","FBLN1","THBS1","LUM","VIM","FN1","COL5A1")
    imm <- marker_score(expr, imm_genes)
    str <- marker_score(expr, str_genes)
    # Standardize and define purity proxy inversely related to stromal+immune
    imm_z <- as.numeric(scale(imm))
    str_z <- as.numeric(scale(str))
    purity_proxy <- - (imm_z + str_z)
    cov$Immune_proxy <- imm
    cov$Stromal_proxy <- str
    cov$Purity_proxy <- purity_proxy
  }

  out_cov <- file.path(opt$outdir, paste0(opt$dataset, ".covariates.tsv"))
  write.table(cov, out_cov, sep="\t", quote=FALSE, row.names=FALSE)
  message("Wrote: ", out_cov)
}

if (identical(environment(), globalenv())){
  tryCatch(main(), error=function(e){
    message("Error: ", e$message)
    quit(status=1)
  })
}
