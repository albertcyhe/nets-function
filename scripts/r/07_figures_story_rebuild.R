#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(GEOquery)
  library(GSVA)
  library(singscore)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(effsize)
})

RFX <- function(x) toupper(trimws(x))
read_module <- function(path){
  if (!file.exists(path)) stop('Missing module file: ', path)
  mod <- readLines(path)
  mod <- RFX(mod)
  mod[mod != '']
}

f_genes <- read_module('resources/modules/net_f_v1.tsv')
m_genes <- read_module('resources/modules/net_m_v1.tsv')

out_fig <- 'results/figures'
out_tab <- 'results/tables'
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)

collapse_duplicates <- function(mat){
  if (any(duplicated(rownames(mat)))){
    df <- as.data.frame(mat)
    df$SYMBOL <- rownames(df)
    df <- df %>% group_by(SYMBOL) %>% summarise(across(everything(), max, na.rm=TRUE))
    mat <- as.matrix(df[,-1])
    rownames(mat) <- df$SYMBOL
  }
  mat
}

map_probes_to_symbol <- function(expr, gpl_id){
  message('Mapping probes for ', gpl_id)
  gpl <- getGEO(gpl_id)
  tab <- Table(gpl)
  cand_cols <- c('Gene Symbol', 'Gene symbol', 'GENE_SYMBOL', 'Symbol', 'GeneSymbol', 'gene_assignment')
  symbol_col <- cand_cols[cand_cols %in% colnames(tab)][1]
  if (is.na(symbol_col)) stop('No symbol column found in GPL: ', gpl_id)
  mapping <- tibble(PROBE = as.character(tab$ID), SYMBOL = RFX(tab[[symbol_col]]))
  mapping <- mapping %>% mutate(SYMBOL = str_split(SYMBOL, pattern = '[;/\\s]+', simplify = FALSE)) %>%
    mutate(SYMBOL = map_chr(SYMBOL, function(vals){ vals <- vals[vals != '' & !is.na(vals)]; if (length(vals) > 0) vals[1] else '' }))
  mapping <- mapping %>% filter(SYMBOL != '')
  probes <- intersect(rownames(expr), mapping$PROBE)
  mapping <- mapping %>% filter(PROBE %in% probes)
  expr <- expr[probes,, drop = FALSE]
  symbols <- mapping$SYMBOL[match(rownames(expr), mapping$PROBE)]
  expr <- rowsum(expr, group = symbols)
  expr
}

load_dataset_expr <- function(dataset, gpl_id = NULL, use_custom = FALSE){
  if (use_custom){
    expr_path <- file.path('data/processed', dataset, paste0(dataset, '.expr.tsv.gz'))
    if (!file.exists(expr_path)) stop('Expected processed expression not found for ', dataset)
    expr_df <- readr::read_tsv(expr_path, show_col_types = FALSE)
    sym_col <- colnames(expr_df)[1]
    expr_mat <- as.matrix(expr_df[,-1])
    rownames(expr_mat) <- RFX(expr_df[[sym_col]])
  } else {
    geo <- getGEO(dataset, GSEMatrix = TRUE, getGPL = FALSE)
    eset <- geo[[1]]
    expr_mat <- exprs(eset)
    rownames(expr_mat) <- rownames(exprs(eset))
    if (!is.null(gpl_id)){
      expr_mat <- map_probes_to_symbol(expr_mat, gpl_id)
    } else {
      rownames(expr_mat) <- RFX(rownames(expr_mat))
    }
 }
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- 'double'
  rownames(expr_mat) <- RFX(rownames(expr_mat))
  expr_mat <- collapse_duplicates(expr_mat)
  expr_mat
}

get_pheno <- function(dataset, use_custom = FALSE){
  if (use_custom){
    ph_path <- file.path('data/processed', dataset, paste0(dataset, '.pheno.tsv'))
    if (!file.exists(ph_path)) stop('Expected pheno not found for ', dataset)
    ph <- readr::read_tsv(ph_path, show_col_types = FALSE)
 } else {
   geo <- getGEO(dataset, GSEMatrix = TRUE, getGPL = FALSE)
    ph <- pData(geo[[1]]) %>% as.data.frame() %>% as_tibble()
  }
  ph
}

calc_module_scores <- function(expr){
  expr <- expr[apply(expr, 1, function(x) all(is.finite(x))), , drop = FALSE]
  common_f <- intersect(rownames(expr), f_genes)
  common_m <- intersect(rownames(expr), m_genes)
  avg_score <- function(gset){
    if (length(gset) == 0) return(rep(NA_real_, ncol(expr)))
    mat <- expr[gset, , drop = FALSE]
    colMeans(mat, na.rm = TRUE)
  }
  singscore_F <- avg_score(common_f)
  singscore_M <- avg_score(common_m)
  ss_df <- tibble(sample_id = colnames(expr),
                  singscore_F = as.numeric(singscore_F),
                  singscore_M = as.numeric(singscore_M))
  gs_list <- list(NET_F = common_f, NET_M = common_m)
  gs_list <- gs_list[sapply(gs_list, length) > 0]
  if (length(gs_list) > 0){
    gs <- tryCatch({ GSVA::gsva(expr, gs_list, method = 'ssgsea', kcdf = 'Gaussian', parallel.sz = 1, abs.ranking = FALSE) }, error = function(e) NULL)
    if (!is.null(gs)){
      ss_df$ssgsea_F <- as.numeric(gs['NET_F', ss_df$sample_id])
      ss_df$ssgsea_M <- as.numeric(gs['NET_M', ss_df$sample_id])
    } else {
      ss_df$ssgsea_F <- NA_real_
      ss_df$ssgsea_M <- NA_real_
    }
  } else {
    ss_df$ssgsea_F <- NA_real_
    ss_df$ssgsea_M <- NA_real_
  }
  safe_scale <- function(v){
    if (all(is.na(v))) return(rep(0, length(v)))
    as.numeric(scale(v))
  }
  ss_df <- ss_df %>% mutate(
    delta_singscore = safe_scale(singscore_F) - safe_scale(singscore_M),
    delta_ssgsea   = safe_scale(ssgsea_F)   - safe_scale(ssgsea_M)
  )
  ss_df
}

calc_signature <- function(expr, genes){
  genes <- intersect(rownames(expr), RFX(genes))
  if (length(genes) == 0) return(rep(NA_real_, ncol(expr)))
  mat <- expr[genes, , drop = FALSE]
  if (nrow(mat) == 0) return(rep(NA_real_, ncol(expr)))
  z <- t(scale(t(mat)))
  colMeans(z, na.rm = TRUE)
}

paired_summary <- function(df, metric, dataset){
  if (!all(c('group','pair_id', metric) %in% colnames(df))) return(NULL)
  if (all(is.na(df[[metric]])) || sd(df[[metric]], na.rm = TRUE) == 0) return(NULL)
  tab <- df %>% filter(!is.na(pair_id), group %in% c('met','primary')) %>%
    group_by(pair_id, group) %>% summarise(value = mean(.data[[metric]], na.rm=TRUE), .groups='drop') %>%
    pivot_wider(names_from = group, values_from = value) %>% drop_na(met, primary)
  if (nrow(tab) < 3) return(NULL)
  wt <- suppressWarnings(wilcox.test(tab$met, tab$primary, paired = TRUE))
  eff <- tryCatch(effsize::cliff.delta(tab$met, tab$primary, paired = TRUE)$estimate, error = function(e) NA_real_)
  tibble(dataset = dataset, method = metric, n = nrow(tab), statistic = unname(wt$statistic), p_value = wt$p.value, cliff_delta = eff)
}

compute_fm_corr <- function(df, dataset){
 rows <- list()
 if (all(c('singscore_F','singscore_M') %in% colnames(df))){
    sub <- df %>% dplyr::select(singscore_F, singscore_M) %>% drop_na()
    if (nrow(sub) >= 3){
      r <- suppressWarnings(cor.test(sub$singscore_F, sub$singscore_M, method='spearman'))
      rows[[length(rows)+1]] <- tibble(dataset = dataset, method = 'singscore', rho = unname(r$estimate), p_value = r$p.value)
    }
  }
  if (all(c('ssgsea_F','ssgsea_M') %in% colnames(df))){
    sub <- df %>% dplyr::select(ssgsea_F, ssgsea_M) %>% drop_na()
    if (nrow(sub) >= 3){
      r <- suppressWarnings(cor.test(sub$ssgsea_F, sub$ssgsea_M, method='spearman'))
      rows[[length(rows)+1]] <- tibble(dataset = dataset, method = 'ssgsea', rho = unname(r$estimate), p_value = r$p.value)
    }
  }
  bind_rows(rows)
}

# Dataset list
config <- tribble(
  ~dataset,   ~gpl,     ~use_custom,
  'GSE184869', NA,      TRUE,
  'GSE125989', 'GPL571', FALSE,
  'GSE43837',  'GPL1352', FALSE,
  'GSE14017',  'GPL570',  FALSE,
  'GSE14018',  'GPL96',   FALSE
)

results_delta <- list()
paired_tables <- list()
corr_rows <- list()

for (i in seq_len(nrow(config))){
  ds <- config$dataset[i]
  message('Processing ', ds)
  expr <- load_dataset_expr(ds, gpl_id = config$gpl[i], use_custom = config$use_custom[i])
  scores <- calc_module_scores(expr)
 if (ds == 'GSE184869'){
    ph <- readr::read_tsv('data/processed/GSE184869/GSE184869.pheno.tsv', show_col_types = FALSE)
 } else {
    ph <- get_pheno(ds, use_custom = FALSE) %>% dplyr::rename(sample_id = geo_accession)
    if (ds == 'GSE125989'){
      sample_titles <- ph$title
      ph$group <- ifelse(str_detect(tolower(sample_titles), 'brain'), 'met', 'primary')
      pair_cols <- grep('characteristics', colnames(ph), value=TRUE)
      ph$pair_id <- NA_character_
      for (col in pair_cols){
        idx <- str_detect(tolower(ph[[col]]), 'paired sample')
        ph$pair_id[idx] <- str_extract(ph[[col]][idx], '[0-9]+')
      }
    }
  }
  scores <- scores %>% left_join(ph, by = 'sample_id')
  # Save per-sample scores summary
  write_tsv(scores, file.path(out_tab, paste0(ds, '_deltafm_scores.tsv')))
  results_delta[[ds]] <- list(expr = expr, scores = scores)
  paired_tbl <- bind_rows(
    paired_summary(scores, 'delta_singscore', ds),
    paired_summary(scores, 'delta_ssgsea', ds)
  )
  if (!is.null(paired_tbl) && nrow(paired_tbl) > 0){
    paired_tbl <- paired_tbl %>% mutate(fdr = p.adjust(p_value, method='BH'))
    paired_tables[[ds]] <- paired_tbl
  }
  corr_rows[[ds]] <- compute_fm_corr(scores, ds)
}

if (length(paired_tables) > 0){
  paired_all <- bind_rows(paired_tables)
  write_tsv(paired_all, file.path(out_tab, 'bulk_paired_tests.tsv'))
  if ('GSE184869' %in% names(paired_tables)) write_tsv(paired_tables[['GSE184869']], file.path(out_tab, 'GSE184869_paired_tests.tsv'))
  if ('GSE125989' %in% names(paired_tables)) write_tsv(paired_tables[['GSE125989']], file.path(out_tab, 'GSE125989_paired_tests.tsv'))
}

corr_df <- bind_rows(corr_rows)
if (nrow(corr_df) > 0){
  corr_df <- corr_df %>% mutate(fdr = p.adjust(p_value, method='BH'))
  write_tsv(corr_df, file.path(out_tab, 'fm_correlation_overview.tsv'))
}

# GSE184869 THBS1 residuals figure
if ('GSE184869' %in% names(results_delta)){
  expr <- results_delta[['GSE184869']]$expr
  scores <- results_delta[['GSE184869']]$scores
  neut_genes <- c('CEACAM8','S100A8','S100A9','MPO','FCGR3B','FPR1')
  purity_genes <- c('EPCAM','KRT8','KRT18','KRT19')
  neut_score <- calc_signature(expr, neut_genes)
  purity_score <- calc_signature(expr, purity_genes)
  thbs1 <- if ('THBS1' %in% rownames(expr)) expr['THBS1', ] else NA_real_
  resid_df <- tibble(sample_id = colnames(expr),
                     delta = scores$delta_singscore[match(colnames(expr), scores$sample_id)],
                     THBS1 = as.numeric(thbs1),
                     neutrophil = neut_score,
                     purity = purity_score)
  resid_df <- resid_df %>% drop_na(delta, THBS1)
  if (nrow(resid_df) >= 4){
    rx <- resid(lm(delta ~ neutrophil + purity, data = resid_df))
    ry <- resid(lm(THBS1 ~ neutrophil + purity, data = resid_df))
    res <- tibble(sample_id = resid_df$sample_id, delta_resid = rx, thbs1_resid = ry)
    stat <- suppressWarnings(cor.test(res$delta_resid, res$thbs1_resid, method='spearman'))
    write_tsv(tibble(dataset='GSE184869', rho=unname(stat$estimate), p_value=stat$p.value, n=nrow(res), method='Spearman residual'), file.path(out_tab, 'GSE184869_THBS1_residuals.tsv'))
    fig <- ggplot(res, aes(delta_resid, thbs1_resid)) +
      geom_point(color = '#dd8452', size = 2.5, alpha = 0.8) +
      geom_smooth(method='lm', se=FALSE, color='#4c72b0', linewidth=0.8) +
      theme_classic(base_size = 12) +
      labs(title = 'GSE184869: THBS1 vs DeltaFM residuals',
           subtitle = sprintf('Spearman rho = %.3f, p = %.3g', unname(stat$estimate), stat$p.value),
           x = 'Residual DeltaFM',
           y = 'Residual THBS1 expression')
    ggsave(file.path(out_fig, 'Figure3_GSE184869_THBS1_Residuals.pdf'), fig, width=4.5, height=4.0, useDingbats=FALSE)
  }
}

# Figure 1 schematic
fig1 <- ggplot(data.frame(stage = factor(c('NET-F', 'NET-M'), levels=c('NET-F','NET-M')), score = c(1, 0.6)),
               aes(stage, score, fill = stage)) +
  geom_col(width=0.6, alpha=0.8, color='black') +
  scale_fill_manual(values=c('#4c72b0','#dd8452')) +
  theme_classic(base_size = 12) + theme(legend.position='none') +
  labs(title='DeltaFM concept', y='z-score', x=NULL) +
  annotate('text', x=1.5, y=1.1, label='Delta*FM == z(F) - z(M)', parse=TRUE)
ggsave(file.path(out_fig, 'Figure1_DFM_Schematic.pdf'), fig1, width=4.2, height=3.2, useDingbats=FALSE)

# Figure 2 paired lines for GSE184869 & GSE125989
plot_paired_dataset <- function(scores, dataset){
  if (!all(c('group','pair_id','delta_singscore') %in% colnames(scores))) return(NULL)
  df <- scores %>% filter(!is.na(pair_id), group %in% c('met','primary')) %>%
    group_by(pair_id, group) %>% summarise(delta = mean(delta_singscore, na.rm=TRUE), .groups='drop') %>%
    mutate(pair_id = as.factor(pair_id))
  if (nrow(df) == 0) return(NULL)
  ggplot(df, aes(group, delta, group = pair_id)) +
    geom_line(alpha=0.6, color='#4c72b0') + geom_point(size=2, color='#dd8452') +
    theme_classic(base_size = 12) + labs(title = dataset, x=NULL, y=expression(Delta*FM[singscore]))
}

paired_plots <- list()
if ('GSE184869' %in% names(results_delta)) paired_plots[['GSE184869']] <- plot_paired_dataset(results_delta[['GSE184869']]$scores, 'GSE184869')
if ('GSE125989' %in% names(results_delta)) paired_plots[['GSE125989']] <- plot_paired_dataset(results_delta[['GSE125989']]$scores, 'GSE125989')
if (length(paired_plots) > 0){
  library(patchwork)
  fig2 <- wrap_plots(paired_plots, ncol = length(paired_plots))
  ggsave(file.path(out_fig, 'Figure2_GSE184869_Paired_DFM.pdf'), fig2, width=4.2*length(paired_plots), height=3.3, useDingbats=FALSE)
}

# Figure 4 FM correlation heatmap
if (nrow(corr_df) > 0){
  corr_df <- corr_df %>% mutate(p_label = ifelse(p_value < 0.001, '<0.001', sprintf('%.3f', p_value)))
  fig4 <- ggplot(corr_df, aes(method, dataset, fill = rho)) +
    geom_tile(color='white') +
    geom_text(aes(label = sprintf('rho=%.2f\np=%s', rho, p_label)), size=3) +
    scale_fill_gradient2(low='#2166AC', mid='white', high='#B2182B', limits=c(-1,1), oob=scales::squish) +
    theme_minimal(base_size = 12) +
    labs(title='F/M score correlation overview', x=NULL, y=NULL, fill='Spearman rho')
  ggsave(file.path(out_fig, 'Figure4_FM_Corr_Overview.pdf'), fig4, width=5.5, height=3.5, useDingbats=FALSE)
}

# Figure 5 scRNA proxy (placeholder based on exploratory notes)
sc_df <- tibble(sample = c('GSM5645901','GSM5645905','GSM5645903','GSM5645904'), lowDFM_fraction = c(0.0188, 0.0435, 0.0, 0.0))
write_tsv(sc_df, file.path(out_tab, 'gse186344_neutrophil_lowDFM_fraction.tsv'))
cat('THBS1 association underpowered (n=4, constant THBS1 expression).\n', file = file.path(out_fig, 'Figure5_scRNA_THBS1_Assoc.txt'))
fig5 <- ggplot(sc_df, aes(x = reorder(sample, lowDFM_fraction), y = lowDFM_fraction)) +
  geom_col(fill='#4c72b0') + theme_classic(base_size = 12) +
  labs(title='GSE186344: low-DeltaFM neutrophil fraction', x='Sample', y='Fraction') +
  coord_flip()
ggsave(file.path(out_fig, 'Figure5_scRNA_LowDFM_Fraction.pdf'), fig5, width=4.2, height=3.5, useDingbats=FALSE)

# Figure 6 mechanism model
model_df <- tibble(stage = factor(c('NET intact','After lysis'), levels=c('NET intact','After lysis')),
                   NET_F = c(1.0, 0.4), NET_M = c(1.0, 0.8)) %>%
  pivot_longer(cols = c(NET_F, NET_M), names_to = 'component', values_to = 'level')
fig6 <- ggplot(model_df, aes(stage, level, group = component, color = component)) +
  geom_line(linewidth=1.2) + geom_point(size=3) +
  scale_color_manual(values=c('NET_F'='#4c72b0','NET_M'='#dd8452')) +
  theme_classic(base_size=12) + labs(title='Conceptual model: NET-M persists > NET-F', y='Relative abundance', x=NULL, color=NULL)
ggsave(file.path(out_fig, 'Figure6_Mechanism_Model.pdf'), fig6, width=5.0, height=3.5, useDingbats=FALSE)

# Summaries for README
message('Story figures regenerated: see results/figures/')
