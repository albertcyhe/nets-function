#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2); library(readr); library(dplyr); library(tidyr)
  library(cowplot); library(ggrepel); library(scales)
})

args <- list(
  assoc = 'results/tables/assoc_summary.tsv',
  serpin = 'results/tables/serpin_scores.tsv',
  footprints = 'results/tables/footprints.tsv',
  thbs1 = 'results/tables/thbs1_cleave_idx.tsv',
  neg = 'results/tables/neg_controls.tsv',
  outdir = 'results/figures'
)

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

assoc <- read_tsv(args$assoc, show_col_types = FALSE) %>%
  mutate(endpoint = factor(endpoint, levels = c('footprint_index','THBS1_log_count','Proteo_DeltaFM'),
                           labels = c('NE/PR3 Footprint','THBS1 (log2 count)','Proteo-ΔFM')),
         stratum = factor(stratum, levels=c('Combined','Brain','NonBrain')))

# Forest plot from Spearman effects with Fisher CI where n>=4
forest_df <- assoc %>% mutate(
  z = atanh(spearman_rho),
  se = 1/sqrt(pmax(n-3, 1)),
  zlo = z - 1.96*se, zhi = z + 1.96*se,
  rlo = tanh(zlo), rhi = tanh(zhi)
)

p_forest <- ggplot(forest_df, aes(x=spearman_rho, y=interaction(stratum, endpoint), color=stratum))+
  geom_point(size=2.2)+
  geom_errorbarh(aes(xmin=rlo, xmax=rhi), height=0.18, size=0.6, alpha=0.8)+
  geom_vline(xintercept = 0, linetype='dashed', color='grey50')+
  scale_x_continuous('Spearman rho (95% CI)', limits=c(-1,1))+
  ylab('Stratum × Endpoint')+
  theme_minimal(base_size = 11)+
  theme(panel.grid.minor = element_blank(), legend.position='bottom')
ggsave(file.path(args$outdir, 'A3_forest_assoc.pdf'), p_forest, width=7.0, height=4.8, useDingbats=FALSE)

# Scatter panels (Serpin vs endpoints)
serpin <- read_tsv(args$serpin, show_col_types = FALSE)
foot <- read_tsv(args$footprints, show_col_types = FALSE)
thbs1 <- read_tsv(args$thbs1, show_col_types = FALSE)

# Merge
neg <- if (file.exists(args$neg)) read_tsv(args$neg, show_col_types = FALSE) else tibble()

scatter <- serpin %>%
  left_join(foot %>% select(dataset, sample, footprint_index), by=c('dataset','sample')) %>%
  left_join(thbs1 %>% select(dataset, sample, Proteo_DeltaFM, THBS1_cleave_idx), by=c('dataset','sample')) %>%
  left_join(neg %>% select(dataset, sample, control_footprint_index, ECM_nc_score), by=c('dataset','sample')) %>%
  # try to bring THBS1 log counts if available (may be NA)
  left_join(read_tsv('data/processed/proteomics/PXD046330/proteo_deltafm.tsv', show_col_types = FALSE) %>% select(Sample, THBS1_log_count) %>% rename(sample=Sample), by='sample') %>%
  mutate(organ = case_when(dataset=='PXD005719' & grepl('BR', toupper(sample)) ~ 'Brain', TRUE ~ 'NonBrain'))

panel_plot <- function(df, ycol, ylab_txt){
  ggplot(df, aes(x=Serpin_score_core, y=.data[[ycol]], color=organ))+
    geom_point(alpha=0.8, size=2)+
    geom_smooth(method='lm', se=FALSE, size=0.7)+
    facet_wrap(~dataset, nrow=1, scales='free_y')+
    labs(x='Serpin score (core, CLR)', y=ylab_txt, color='Stratum')+
    theme_classic(base_size = 11)+
    theme(legend.position='bottom')
}

p1 <- panel_plot(scatter, 'footprint_index', 'NE/PR3 Footprint index')
p2 <- panel_plot(scatter, 'THBS1_log_count', 'THBS1 (log2 count)')
p3 <- panel_plot(scatter, 'Proteo_DeltaFM', 'Proteo-ΔFM')
p4 <- panel_plot(scatter, 'THBS1_cleave_idx', 'THBS1 cleavage index')
pc1 <- panel_plot(scatter, 'control_footprint_index', 'Control footprint index (acidic P1)')
pc2 <- panel_plot(scatter, 'ECM_nc_score', 'ECM negative-control score (CLR)')

ggsave(file.path(args$outdir, 'A3_scatter_serpin_vs_footprint.pdf'), p1, width=8.5, height=2.8, useDingbats=FALSE)
ggsave(file.path(args$outdir, 'A3_scatter_serpin_vs_THBS1.pdf'), p2, width=8.5, height=2.8, useDingbats=FALSE)
ggsave(file.path(args$outdir, 'A3_scatter_serpin_vs_deltaFM.pdf'), p3, width=8.5, height=2.8, useDingbats=FALSE)
ggsave(file.path(args$outdir, 'A3_scatter_serpin_vs_thbs1_cleave.pdf'), p4, width=8.5, height=2.8, useDingbats=FALSE)
ggsave(file.path(args$outdir, 'A3_scatter_serpin_vs_ctrl_footprint.pdf'), pc1, width=8.5, height=2.8, useDingbats=FALSE)
ggsave(file.path(args$outdir, 'A3_scatter_serpin_vs_ecm_nc.pdf'), pc2, width=8.5, height=2.8, useDingbats=FALSE)

message('Wrote figures to ', args$outdir)
