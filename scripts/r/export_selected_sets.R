#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(msigdbr)
})

main <- function(){
  dir.create('resources/pathways', recursive = TRUE, showWarnings = FALSE)
  # Load cached if available, else query
  msig_c2r <- if (file.exists('resources/pathways/msig_c2_reactome.rds')) readRDS('resources/pathways/msig_c2_reactome.rds') else msigdbr(species='Homo sapiens', collection='C2', subcollection='REACTOME')
  msig_c5  <- if (file.exists('resources/pathways/msig_c5_gobp.rds')) readRDS('resources/pathways/msig_c5_gobp.rds') else msigdbr(species='Homo sapiens', collection='C5', subcollection='GO:BP')

  re_names <- unique(msig_c2r$gs_name)
  go_names <- unique(msig_c5$gs_name)
  re_sel <- re_names[grepl('P38|P-38', re_names, ignore.case = TRUE)]
  go_sel <- go_names[grepl('INTERLEUKIN_1|IL1', go_names, ignore.case = TRUE)]

  tab <- rbind(
    data.frame(source='REACTOME', set=re_sel, stringsAsFactors = FALSE),
    data.frame(source='GO:BP', set=go_sel, stringsAsFactors = FALSE)
  )
  out <- 'resources/pathways/selected_sets.tsv'
  write.table(tab, out, sep='\t', quote=FALSE, row.names=FALSE)
  message('Wrote: ', out)
}

if (identical(environment(), globalenv())){
  tryCatch(main(), error=function(e){ message('Error: ', e$message); quit(status=1) })
}

