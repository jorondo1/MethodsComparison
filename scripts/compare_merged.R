# Compare the phyloseq objects built from the merged tables (scripts/utils/parse_merged.R, merge_tables.sh outputs)
# with those built from the per-sample files, the way 0_proc_data.R did before the merged tables
# (community_functions.R pinned at misc_scripts b86aac6, so the reference can't move).
# For every dataset x database: samples, taxa, tax_table, otu_table, sample sums and alpha diversity of the raw objects.
#
# Usage (from $MC): Rscript scripts/compare_merged.R [dataset ...]   (default: every dataset with sample data, + Maple)
# On the HPC: module load StdEnv/2023 r/4.4.0
# Writes Out/merged_tables/comparison.tsv and prints a summary.

suppressMessages({
  library(phyloseq); library(mgx.tools); library(tidyverse); library(magrittr); library(vegan)
})

REF_FUNCTIONS <- 'https://raw.githubusercontent.com/jorondo1/misc_scripts/b86aac6/community_functions.R'
ref <- new.env()
ref_file <- tempfile(fileext = '.R'); download.file(REF_FUNCTIONS, ref_file, quiet = TRUE)
invisible(capture.output(sys.source(ref_file, envir = ref)))
new <- new.env()
sys.source('scripts/utils/parse_merged.R', envir = new)
source('scripts/utils/sample_data.R')

# Reference: the per-database blocks of the former meta_parsing() (0_proc_data.R before the merged tables), one at a time
ref_parse <- function(dsName, db) {
  assemble <- function(tbl, ...) assemble_phyloseq(tbl, samData, seqdepth_plot = FALSE, onlySpecies = TRUE, ...)
  path <- function(...) file.path('data', dsName, ...)
  parse_MPA <- ref$parse_MPA; parse_SM <- ref$parse_SM; species_glom <- ref$species_glom
  parse_GTDB_lineages <- ref$parse_GTDB_lineages; parse_genbank_lineages <- ref$parse_genbank_lineages
  {
    if (db == 'MOTUS3') {
      parse_MPA(MPA_files = path("MOTUS/*_profile.txt"), column_names = c('mOTU', 'Taxonomy', 'NCBI', 'Abundance'),
                mOTUs_data = TRUE) %>% assemble()
    } else if (db == 'MOTUS4') {
      parse_MPA(MPA_files = path("MOTUS4/*_profile.txt"), column_names = c('mOTU', 'Taxonomy', 'Abundance'),
                mOTUs_data = TRUE) %>% assemble()
    } else if (str_detect(db, '^SM_')) {
      lineage_path <- path(db, paste0(db, '_lineages.csv'))
      left_join(
        parse_SM(path(db, '*_gather.csv')),
        if (str_detect(db, 'gtdb')) parse_GTDB_lineages(lineage_path) else parse_genbank_lineages(lineage_path),
        by = 'genome') %>%
        species_glom() %>% assemble()
    } else if (str_detect(db, '^KB')) {
      parse_MPA(MPA_files = path(db, '*/*_bracken/*_bracken_S.MPA.TXT')) %>% assemble()
    } else if (str_detect(db, '^MPA_')) {
      parse_MPA(MPA_files = path(db, '*/*_profile.txt'), column_names = c('Taxonomy', 'NCBI', 'Abundance', 'Void'),
                MPA_data = TRUE) %>%
        assemble(min_sample_size = 0, min_taxa_count = 0)
    } else if (str_detect(db, '^SYLPH_')) {
      # same reader on both sides (sylph was already merged): checks the shared code path only
      new$parse_sylph(path(db)) %>% species_glom() %>% assemble()
    } else stop('unknown database ', db)
  }
}

otu <- function(ps) { m <- as(otu_table(ps), 'matrix'); if (!taxa_are_rows(ps)) m <- t(m); m }

compare_ps <- function(a, b) { # a = reference, b = merged
  if (is.null(a) || is.null(b)) return(tibble(status = 'missing object'))
  ma <- otu(a); mb <- otu(b)
  s_common <- intersect(colnames(ma), colnames(mb)); t_common <- intersect(rownames(ma), rownames(mb))
  ta <- as(tax_table(a), 'matrix'); tb <- as(tax_table(b), 'matrix')
  da <- ma[t_common, s_common, drop = FALSE]; db <- mb[t_common, s_common, drop = FALSE]
  # values absent from one side's taxa count as differences
  sums_a <- colSums(ma)[s_common]; sums_b <- colSums(mb)[s_common]
  rel <- function(x, y) ifelse(pmax(abs(x), abs(y)) > 0, abs(x - y) / pmax(abs(x), abs(y)), 0)
  shannon_a <- vegan::diversity(t(ma[, s_common, drop = FALSE])); shannon_b <- vegan::diversity(t(mb[, s_common, drop = FALSE]))
  rich_a <- colSums(ma[, s_common, drop = FALSE] > 0); rich_b <- colSums(mb[, s_common, drop = FALSE] > 0)
  tibble(
    n_samples_ref = ncol(ma), n_samples_new = ncol(mb),
    samples_only_ref = length(setdiff(colnames(ma), colnames(mb))),
    samples_only_new = length(setdiff(colnames(mb), colnames(ma))),
    n_taxa_ref = nrow(ma), n_taxa_new = nrow(mb),
    taxa_only_ref = length(setdiff(rownames(ma), rownames(mb))),
    taxa_only_new = length(setdiff(rownames(mb), rownames(ma))),
    tax_table_equal = identical(colnames(ta), colnames(tb)) &&
      identical(unname(ta[t_common, , drop = FALSE]), unname(tb[t_common, , drop = FALSE])),
    otu_max_abs_diff = if (length(da)) max(abs(da - db)) else NA,
    otu_max_rel_diff = if (length(da)) max(rel(da, db)) else NA,
    sum_ratio_ref_new_median = median(sums_a / sums_b),
    sample_sums_max_rel_diff = max(rel(sums_a, sums_b)),
    shannon_max_abs_diff = max(abs(shannon_a - shannon_b)),
    richness_max_abs_diff = max(abs(rich_a - rich_b))
  ) %>%
    mutate(status = case_when(
      samples_only_ref + samples_only_new + taxa_only_ref + taxa_only_new > 0 | !tax_table_equal ~ 'DIFFERENT',
      otu_max_rel_diff == 0 ~ 'identical',
      otu_max_rel_diff < 1e-9 ~ 'equal (float)',
      TRUE ~ 'DIFFERENT'))
}

sam_data_all <- load_sample_data()
datasets <- commandArgs(trailingOnly = TRUE)
if (length(datasets) == 0) datasets <- c('P19_Saliva', 'P19_Gut', 'Moss', 'NAFLD', 'AD_Skin', 'RA_Gut', 'Bee', 'PD', 'Maple')

# Maple isn't in 0_proc_data.R yet: compare with every sample of its sample tsv
sample_tsv <- function(ds) {
  l <- readLines('scripts/src/datasets.tsv') %>% grep('^#', ., invert = TRUE, value = TRUE) %>% strsplit('\t')
  p <- Filter(function(x) x[1] == ds, l)[[1]][2] %>% str_replace('\\$ILAFORES', Sys.getenv('ILAFORES', '/net/nfs-bio/jbod2/def-ilafores'))
  read_tsv(p, col_names = FALSE, show_col_types = FALSE, col_types = cols(.default = 'c'))$X1
}

results <- list()
for (ds in datasets) {
  samData <- if (!is.null(sam_data_all[[ds]])) sam_data_all[[ds]] else data.frame(row.names = sample_tsv(ds), all = 1)
  message('### ', ds)
  new_ps <- new$meta_parsing(ds, samData)
  dirs <- new$tool_dirs(ds)
  dbs <- c(if (length(dirs$motus3)) 'MOTUS3', if (length(dirs$motus4)) 'MOTUS4',
           dirs$sourmash, dirs$sylph[file.exists(file.path('data', ds, dirs$sylph, 'counts.tsv'))],
           dirs$kraken, dirs$metaphlan)
  for (db in dbs) {
    message('  ', db)
    t0 <- Sys.time()
    ref_ps <- tryCatch(ref_parse(ds, db), error = function(e) { message('    reference failed: ', conditionMessage(e)); NULL })
    results[[paste(ds, db)]] <- compare_ps(ref_ps, new_ps[[db]]) %>%
      mutate(dataset = ds, database = db, .before = 1)
    message(sprintf('    %s (%.0f s)', results[[paste(ds, db)]]$status, as.numeric(difftime(Sys.time(), t0, units = 'secs'))))
  }
}

res <- bind_rows(results)
dir.create('Out/merged_tables', recursive = TRUE, showWarnings = FALSE)
write_tsv(res, 'Out/merged_tables/comparison.tsv')
cat('\n'); print(count(res, status))
cat('\nNot identical:\n')
res %>% filter(status != 'identical') %>%
  select(dataset, database, status, samples_only_ref, samples_only_new, taxa_only_ref, taxa_only_new,
         tax_table_equal, otu_max_rel_diff, sum_ratio_ref_new_median) %>%
  as.data.frame() %>% print()
cat('\nWritten to Out/merged_tables/comparison.tsv\n')
