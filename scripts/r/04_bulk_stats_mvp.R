#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(effsize)
})

read_tsv_dt <- function(path){
  fread(path, sep='\t', header=TRUE, data.table=FALSE, check.names=FALSE)
}

detect_groups_pairs_gse125989 <- function(pheno){
  df <- pheno
  df$group <- NA_character_
  if ('title' %in% colnames(df)){
    tl <- tolower(df$title)
    df$group[grepl('brain metast', tl)] <- 'met'
    df$group[grepl('primary', tl) & is.na(df$group)] <- 'primary'
  }
  # paired sample id from characteristics columns
  pair <- apply(df, 1, function(r){
    vals <- as.character(r)
    m <- grep('paired sample:', tolower(vals), value=TRUE)
    if (length(m)>0){
      gsub('.*paired sample:\\s*', '', tolower(m[1]))
    } else NA_character_
  })
  df$pair_id <- pair
  df
}

detect_groups_pairs_gse184869 <- function(pheno){
  df <- pheno
  sid <- df$sample_id
  grp <- rep(NA_character_, length(sid))
  pid <- rep(NA_character_, length(sid))
  # Pattern 1: nM_RCS / nP_RCS
  m <- regexec('^([0-9]+)([MP])_RCS$', sid)
  mm <- regmatches(sid, m)
  for (i in seq_along(mm)){
    if (length(mm[[i]])==3){
      pid[i] <- mm[[i]][2]
      grp[i] <- ifelse(mm[[i]][3]=='M','met','primary')
    }
  }
  # Pattern 2: BMxx / BPxx
  idx <- which(is.na(grp))
  m2 <- regexec('^B([MP])(\\d+)', sid[idx])
  mm2 <- regmatches(sid[idx], m2)
  for (k in seq_along(idx)){
    j <- idx[k]
    if (length(mm2[[k]])==3){
      pid[j] <- mm2[[k]][3]
      grp[j] <- ifelse(mm2[[k]][2]=='M','met','primary')
    }
  }
  df$group <- grp
  df$pair_id <- pid
  df
}

partial_cor_resid <- function(x, y, covars){
  # residualize x and y on covars, then Pearson/Spearman on residuals
  dat <- data.frame(x=x, y=y, covars)
  dat <- dat[complete.cases(dat),,drop=FALSE]
  if (nrow(dat) < 3) return(list(r=NA, p=NA, n=nrow(dat)))
  rx <- resid(lm(x ~ ., data=dat))
  ry <- resid(lm(y ~ ., data=dat))
  ct <- suppressWarnings(cor.test(rx, ry, method='spearman'))
  list(r=unname(ct$estimate), p=ct$p.value, n=length(rx))
}

collapse_duplicated_rows <- function(expr){
  if (any(duplicated(rownames(expr)))){
    expr <- t(apply(expr, 2, function(col){ tapply(col, INDEX = rownames(expr), FUN = function(v) max(v, na.rm=TRUE)) }))
    expr <- t(expr)
  }
  expr
}

paired_test <- function(x1, x2){
  keep <- is.finite(x1) & is.finite(x2)
  x1 <- x1[keep]; x2 <- x2[keep]
  if (length(x1) < 3) return(list(stat=NA, p=NA, n=length(x1)))
  wt <- suppressWarnings(wilcox.test(x2, x1, paired=TRUE, alternative='two.sided'))
  list(stat=unname(wt$statistic), p=wt$p.value, n=length(x1))
}

pick_cov <- function(tab, name, fallback=NULL){
  # prefer exact, then plural, then proxy
  for (cand in c(name, paste0(name,'s'), paste0(name,'_proxy'), fallback)){
    if (!is.null(cand) && cand %in% colnames(tab)) return(tab[[cand]])
  }
  rep(NA_real_, nrow(tab))
}

infer_covar_sources <- function(ds, tab, outdir_base='data/processed'){
  # Neutrophil source
  neut_src <- if ('Neutrophil' %in% colnames(tab)) 'MCPcounter:Neutrophil'
              else if ('Neutrophils' %in% colnames(tab)) 'MCPcounter:Neutrophils'
              else if ('Neutrophil_proxy' %in% colnames(tab)) 'Proxy:marker_score'
              else 'Unknown'
  # Purity source
  purity_src <- 'Proxy:marker_score'
  if ('TumorPurity' %in% colnames(tab)){
    # Inspect platform file to distinguish Affy vs RNA-seq derived
    plat_file <- file.path(outdir_base, ds, paste0(ds, '.platform.tsv'))
    if (file.exists(plat_file)){
      pf <- tryCatch(read.delim(plat_file, check.names = FALSE), error=function(e) NULL)
      if (!is.null(pf)){
        if ('Series_platform_id' %in% colnames(pf)) purity_src <- 'tidyestimate:Affymetrix'
        if ('platform' %in% colnames(pf) && grepl('rna', tolower(pf$platform[1]))) purity_src <- 'tidyestimate:ESTIMATEScore-derived'
      } else purity_src <- 'tidyestimate'
    } else purity_src <- 'tidyestimate'
  }
  list(neutrophil_source=neut_src, purity_source=purity_src)
}

main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop('Usage: Rscript 04_bulk_stats_mvp.R <dataset>')
  ds <- args[[1]]

  expr_path <- file.path('data/processed', ds, paste0(ds, '.expr.tsv.gz'))
  pheno_path <- file.path('data/processed', ds, paste0(ds, '.pheno.tsv'))
  scores_path <- file.path('data/processed', ds, paste0(ds, '.scores.tsv'))
  cov_path <- file.path('data/processed', ds, paste0(ds, '.covariates.tsv'))
  outdir <- file.path('results', 'tables')
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  expr <- read_tsv_dt(expr_path)
  rownames(expr) <- expr[[1]]; expr <- expr[,-1,drop=FALSE]
  # ensure unique rownames
  expr <- collapse_duplicated_rows(as.matrix(expr))
  pheno <- read_tsv_dt(pheno_path)
  scores <- read_tsv_dt(scores_path)
  cov <- read_tsv_dt(cov_path)

  # annotate groups/pairs (dataset specific logic)
  anno <- pheno
  if (ds == 'GSE125989') anno <- detect_groups_pairs_gse125989(pheno)
  if (ds == 'GSE184869') anno <- detect_groups_pairs_gse184869(pheno)

  # Merge scores + covariates + anno
  tab <- merge(scores, cov, by='sample_id', all.x=TRUE)
  cols <- intersect(c('sample_id','group','pair_id'), colnames(anno))
  tab <- merge(tab, anno[, cols, drop=FALSE], by='sample_id', all.x=TRUE)
  # Recompute ΔFM within-dataset to avoid stored column issues
  zfun <- function(x){ x <- as.numeric(x); (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE) }
  if (all(c('singscore_F','singscore_M') %in% colnames(tab))){
    tab$deltaFM_singscore_re <- zfun(tab$singscore_F) - zfun(tab$singscore_M)
  }
  if (all(c('ssgsea_F','ssgsea_M') %in% colnames(tab))){
    tab$deltaFM_ssgsea_re <- zfun(tab$ssgsea_F) - zfun(tab$ssgsea_M)
  }

  # paired analysis (primary vs met within pair)
  paired_res <- data.frame()
  if ('pair_id' %in% colnames(tab) && 'group' %in% colnames(tab)){
    # For singscore (aggregate duplicates by mean)
    metric <- if ('deltaFM_singscore_re' %in% colnames(tab)) 'deltaFM_singscore_re' else 'deltaFM_singscore'
    agg <- aggregate(tab[[metric]] ~ tab$pair_id + tab$group, FUN=function(x) mean(x, na.rm=TRUE))
    colnames(agg) <- c('pair_id','group','delta')
    wide <- reshape(agg, idvar='pair_id', timevar='group', direction='wide')
    d1 <- wide$`delta.met` - wide$`delta.primary`
    d1 <- d1[is.finite(d1)]
    wt1 <- suppressWarnings(wilcox.test(d1, mu=0))
    es <- tryCatch({ effsize::cliff.delta(d1, rep(0, length(d1)))$estimate }, error=function(e) NA)
    src <- infer_covar_sources(ds, tab)
    paired_res <- rbind(paired_res, data.frame(dataset=ds, method='singscore', stat=unname(wt1$statistic), p=wt1$p.value, n=length(d1), effect=es, neutrophil_source=src$neutrophil_source, purity_source=src$purity_source))
    # For ssGSEA if available
    if (('deltaFM_ssgsea_re' %in% colnames(tab)) || ('deltaFM_ssgsea' %in% colnames(tab))){
      metric2 <- if ('deltaFM_ssgsea_re' %in% colnames(tab)) 'deltaFM_ssgsea_re' else 'deltaFM_ssgsea'
      agg2 <- aggregate(tab[[metric2]] ~ tab$pair_id + tab$group, FUN=function(x) mean(x, na.rm=TRUE))
      colnames(agg2) <- c('pair_id','group','delta')
      wide2 <- reshape(agg2, idvar='pair_id', timevar='group', direction='wide')
      d2 <- wide2$`delta.met` - wide2$`delta.primary`
      d2 <- d2[is.finite(d2)]
      wt2 <- suppressWarnings(wilcox.test(d2, mu=0))
      es2 <- tryCatch({ effsize::cliff.delta(d2, rep(0, length(d2)))$estimate }, error=function(e) NA)
      src <- infer_covar_sources(ds, tab)
      paired_res <- rbind(paired_res, data.frame(dataset=ds, method='ssgsea', stat=unname(wt2$statistic), p=wt2$p.value, n=length(d2), effect=es2, neutrophil_source=src$neutrophil_source, purity_source=src$purity_source))
    }
  }
  if (nrow(paired_res)>0){ paired_res$fdr <- p.adjust(paired_res$p, method='BH') }
  write.table(paired_res, file.path(outdir, paste0(ds, '_paired_tests.tsv')), sep='\t', quote=FALSE, row.names=FALSE)

  # THBS1 association (partial correlation controlling neutrophil/purity proxies)
  # Map THBS1 if needed (Ensembl case)
  if (!('THBS1' %in% rownames(expr))){
    suppressWarnings({
      if (requireNamespace('AnnotationDbi', quietly=TRUE) && requireNamespace('org.Hs.eg.db', quietly=TRUE)){
        rn <- rownames(expr)
        if (grepl('^ENSG', rn[1])){
          ens <- gsub('\\..*$', '', rn)
          map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=unique(ens), columns=c('SYMBOL'), keytype='ENSEMBL')
          map <- map[!is.na(map$SYMBOL) & map$SYMBOL!='', c('ENSEMBL','SYMBOL')]
          idx <- match(ens, map$ENSEMBL)
          sy <- map$SYMBOL[idx]
          sy[is.na(sy)] <- rn[is.na(sy)]
          rownames(expr) <- toupper(sy)
          expr <- collapse_duplicated_rows(expr)
        }
      }
    })
  }
  thbs1 <- if ('THBS1' %in% rownames(expr)) as.numeric(expr['THBS1', tab$sample_id]) else rep(NA_real_, nrow(tab))

  covars <- data.frame(Neutrophil = pick_cov(tab, 'Neutrophil'),
                       Purity = pick_cov(tab, 'TumorPurity', fallback='Purity_proxy'))
  metric <- if ('deltaFM_singscore_re' %in% colnames(tab)) tab$deltaFM_singscore_re else tab$deltaFM_singscore
  pc1 <- partial_cor_resid(metric, thbs1, covars)
  src <- infer_covar_sources(ds, tab)
  res_cor <- data.frame(dataset=ds, method='singscore', target='THBS1', r=pc1$r, p=pc1$p, n=pc1$n, neutrophil_source=src$neutrophil_source, purity_source=src$purity_source)
  if ('deltaFM_ssgsea_re' %in% colnames(tab) || 'deltaFM_ssgsea' %in% colnames(tab)){
    metric2 <- if ('deltaFM_ssgsea_re' %in% colnames(tab)) tab$deltaFM_ssgsea_re else tab$deltaFM_ssgsea
    pc2 <- partial_cor_resid(metric2, thbs1, covars)
    res_cor <- rbind(res_cor, data.frame(dataset=ds, method='ssgsea', target='THBS1', r=pc2$r, p=pc2$p, n=pc2$n, neutrophil_source=src$neutrophil_source, purity_source=src$purity_source))
  }
  if (nrow(res_cor)>0){ res_cor$fdr <- p.adjust(res_cor$p, method='BH') }
  write.table(res_cor, file.path(outdir, paste0(ds, '_thbs1_partial_cor.tsv')), sep='\t', quote=FALSE, row.names=FALSE)
}

if (identical(environment(), globalenv())){
  tryCatch(main(), error=function(e){
    message('Error: ', e$message)
    quit(status=1)
  })
}
