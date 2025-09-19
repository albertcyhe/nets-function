#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

outdir_fig <- "results/figures"
dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)

read_scores <- function(ds){
  f <- file.path("data/processed", ds, paste0(ds, ".scores.tsv"))
  if (!file.exists(f)) stop("scores not found ", f)
  fread(f)
}
read_pheno <- function(ds){
  f <- file.path("data/processed", ds, paste0(ds, ".pheno.tsv"))
  fread(f)
}

zs <- function(x){ (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE) }

detect_pairs_gse125989 <- function(ph){
  ph$group <- ifelse(grepl("brain metast", tolower(ph$title)), "met", "primary")
  pid <- apply(ph,1,function(r){ m<-grep("paired sample:",tolower(as.character(r)),value=TRUE); if(length(m)>0) sub(".*paired sample: ","",m[1]) else NA })
  ph$pair_id <- pid
  ph
}

detect_pairs_gse184869 <- function(ph){
  sid <- ph$sample_id
  grp <- rep(NA_character_, length(sid)); pid <- rep(NA_character_, length(sid))
  m <- regexec('^([0-9]+)([MP])_RCS$', sid); mm <- regmatches(sid,m)
  for(i in seq_along(mm)) if(length(mm[[i]])==3){ pid[i]<-mm[[i]][2]; grp[i]<-ifelse(mm[[i]][3]=='M','met','primary') }
  idx <- which(is.na(grp))
  m2 <- regexec('^B([MP])(\\d+)', sid[idx]); mm2 <- regmatches(sid[idx], m2)
  for(k in seq_along(idx)) if(length(mm2[[k]])==3){ j<-idx[k]; pid[j]<-mm2[[k]][3]; grp[j]<-ifelse(mm2[[k]][2]=='M','met','primary') }
  ph$group <- grp; ph$pair_id <- pid; ph
}

plot_paired_lines <- function(ds){
  sc <- read_scores(ds)
  ph <- read_pheno(ds)
  if(ds=="GSE125989") ph <- detect_pairs_gse125989(ph)
  if(ds=="GSE184869") ph <- detect_pairs_gse184869(ph)
  tab <- merge(sc, ph[,c("sample_id","group","pair_id")], by="sample_id", all.x=TRUE)
  # recompute ΔFM
  tab$deltaFM_singscore_re <- zs(tab$singscore_F) - zs(tab$singscore_M)
  keep <- !is.na(tab$pair_id) & !is.na(tab$group)
  tab <- tab[keep,]
  if(nrow(tab)==0) return(NULL)
  # long to plot
  agg <- tab %>% dplyr::select(pair_id, group, deltaFM_singscore_re)
  wide <- reshape(agg, idvar='pair_id', timevar='group', direction='wide')
  long <- tidyr::pivot_longer(as.data.frame(wide), cols = c(`deltaFM_singscore_re.met`,`deltaFM_singscore_re.primary`),
                               names_to="grp", values_to="dfm")
  long$group <- ifelse(grepl("met$", long$grp), "met", "primary")
  long$pair_id <- rep(wide$pair_id, each=2)
  p <- ggplot(long, aes(x=group, y=dfm, group=pair_id))+
    geom_line(alpha=0.6)+geom_point(size=2)+
    theme_classic()+labs(title=paste0(ds, " ΔFM (singscore) paired"), x="", y=expression(Delta*FM))
  ggsave(file.path(outdir_fig, paste0("Figure2_", ds, "_Paired_DFM.pdf")), p, width=4.5, height=4)
}

plot_thbs1_scatter_gse184869 <- function(){
  # residual scatter plot using recomputed ΔFM and THBS1 expression (log scale), controlling covariates -> partial via residuals
  ds <- "GSE184869"
  sc <- as.data.frame(read_scores(ds)); ph <- as.data.frame(read_pheno(ds))
  exprf <- file.path("data/processed", ds, paste0(ds, ".expr.tsv.gz"))
  ex <- fread(exprf)
  rn <- ex[[1]]; ex <- as.data.frame(ex[,-1]); rownames(ex) <- rn
  # coerce to numeric and collapse duplicate Ensembl rows (pre-mapping)
  exm <- as.matrix(sapply(ex, function(x) suppressWarnings(as.numeric(x))))
  rownames(exm) <- rownames(ex)
  if(any(duplicated(rownames(exm)))){
    idx_list <- split(seq_len(nrow(exm)), rownames(exm))
    exm <- do.call(rbind, lapply(idx_list, function(ix){
      apply(exm[ix,,drop=FALSE], 2, function(v) max(v, na.rm=TRUE))
    }))
  }
  ex <- as.data.frame(exm)
  # map Ensembl -> SYMBOL if needed
  if(grepl('^ENSG', rownames(ex)[1])){
    suppressMessages({library(AnnotationDbi); library(org.Hs.eg.db)})
    ens <- gsub("\\..*$","", rownames(ex))
    map <- AnnotationDbi::select(org.Hs.eg.db, keys=unique(ens), columns=c("SYMBOL"), keytype="ENSEMBL")
    sy <- map$SYMBOL[match(ens, map$ENSEMBL)]; sy[is.na(sy)] <- rownames(ex)
    sym <- toupper(sy)
    # build a new matrix collapsed by SYMBOL (col-wise max)
    exm <- as.matrix(sapply(ex, function(x) suppressWarnings(as.numeric(x))))
    rownames(exm) <- sym
    idx_list <- split(seq_len(nrow(exm)), rownames(exm))
    exm2 <- do.call(rbind, lapply(idx_list, function(ix){
      apply(exm[ix,,drop=FALSE], 2, function(v) max(v, na.rm=TRUE))
    }))
    ex <- as.data.frame(exm2)
    }
  if(!("THBS1" %in% rownames(ex))) return(NULL)
  df <- sc %>% dplyr::mutate(dfm = zs(singscore_F) - zs(singscore_M)) %>% dplyr::select(sample_id, dfm)
  covf <- file.path("data/processed", ds, paste0(ds, ".covariates.tsv"))
  cov <- as.data.frame(fread(covf))
  cov <- cov %>% dplyr::select(sample_id, Neutrophil = matches("Neutrophil|Neutrophils"), TumorPurity = matches("TumorPurity|Purity_proxy"))
  dd <- df %>% left_join(cov, by="sample_id")
  # residualize
  fit_x <- lm(dfm ~ Neutrophil + TumorPurity, data=dd)
  rx <- resid(fit_x)
  y <- as.numeric(ex["THBS1", dd$sample_id])
  fit_y <- lm(y ~ dd$Neutrophil + dd$TumorPurity)
  ry <- resid(fit_y)
  d2 <- data.frame(rx=rx, ry=ry)
  r <- suppressWarnings(cor.test(d2$rx, d2$ry, method='spearman'))
  p <- ggplot(d2, aes(rx, ry)) + geom_point(alpha=0.7) + geom_smooth(method='lm', se=FALSE, color='red') +
    theme_classic() + labs(title=paste0("GSE184869 THBS1 vs DeltaFM (residuals)\n",
                                         sprintf("rho=%.3f, p=%.3g", unname(r$estimate), r$p.value)),
                           x="DeltaFM residual (DeltaFM ~ N + Purity)", y="THBS1 residual (THBS1 ~ N + Purity)")
  ggsave(file.path(outdir_fig, "Figure3_GSE184869_THBS1_Residuals.pdf"), p, width=4.8, height=4, useDingbats=FALSE)
}

plot_fm_corr_overview <- function(){
  f <- "results/tables/fm_correlation_overview.tsv"
  if(!file.exists(f)) return(NULL)
  d <- fread(f)
  d$label <- sprintf("rho=%.2f\nFDR=%.3g", d$rho, d$fdr)
  p <- ggplot(d, aes(x=method, y=dataset, fill=rho)) +
    geom_tile(color='white') + geom_text(aes(label=label), size=3) + scale_fill_gradient2(low='blue', high='red', mid='white') +
    theme_minimal() + labs(title="F/M correlation overview (Spearman)")
  ggsave(file.path(outdir_fig, "Figure4_FM_Corr_Overview.pdf"), p, width=6, height=3.5)
}

make_schematic <- function(){
  # simple ΔFM concept schematic using ggplot
  df <- data.frame(x=c(1,2), y=c(1,1), lab=c("NET-F (function)", "NET-M (marker)"))
  p <- ggplot(df, aes(x,y,label=lab, fill=lab)) + geom_tile(width=0.8, height=0.5, alpha=0.6) +
    geom_segment(aes(x=1.4,xend=1.6,y=1.2,yend=1.2), arrow=arrow(length=unit(0.15,'cm')))+
    annotate("text", x=1.5, y=1.35, label="DeltaFM = z(F) - z(M)")+
    theme_void() + theme(legend.position='none') + xlim(0.5,2.5)+ylim(0.5,1.8)
  ggsave(file.path(outdir_fig, "Figure1_DFM_Schematic.pdf"), p, width=5, height=3, useDingbats=FALSE)
}

make_scRNA_summary <- function(){
  fracf <- "results/tables/gse186344_neutrophil_lowDFM_fraction.tsv"
  assocf <- "results/tables/gse186344_subpop_thbs1_assoc.tsv"
  if(!file.exists(fracf)) return(NULL)
  fr <- fread(fracf)
  p1 <- ggplot(fr, aes(x=sample, y=lowDFM_neut_fraction)) + geom_col() + theme_minimal() +
    labs(title="Low-ΔFM neutrophil fraction per sample", y="fraction", x="sample") + coord_flip()
  ggsave(file.path(outdir_fig, "Figure5_scRNA_LowDFM_Fraction.pdf"), p1, width=5, height=3.5)
  if(file.exists(assocf)){
    as <- fread(assocf)
    con <- ifelse(is.na(as$rho[1]), "THBS1 association NA (constant or missing)", sprintf("rho=%.3f, p=%.3g, n=%d", as$rho[1], as$p[1], as$n[1]))
    cat(con, file = file.path(outdir_fig, "Figure5_scRNA_THBS1_Assoc.txt"))
  }
}

make_model_figure <- function(){
  # conceptual model: after NETs breakdown, marker (M) persists more than function (F)
  df <- data.frame(stage=c("NET intact","After lysis"), F=c(1.0, 0.4), M=c(1.0, 0.8))
  dfl <- tidyr::pivot_longer(df, cols=c(F,M), names_to="comp", values_to="level")
  p <- ggplot(dfl, aes(stage, level, group=comp, color=comp)) + geom_line(size=1.2) + geom_point(size=2) +
    theme_classic() + labs(title="Proposed model: NET-M persists > NET-F in microenvironment", y="relative level")
  ggsave(file.path(outdir_fig, "Figure6_Mechanism_Model.pdf"), p, width=5.5, height=3.5)
}

main <- function(){
  make_schematic()
  suppressWarnings(plot_paired_lines("GSE184869"))
  suppressWarnings(plot_paired_lines("GSE125989"))
  suppressWarnings(plot_thbs1_scatter_gse184869())
  plot_fm_corr_overview()
  make_scRNA_summary()
  make_model_figure()
  message("Figures written under ", outdir_fig)
}

if (identical(environment(), globalenv())) main()
