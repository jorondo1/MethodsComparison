# Master script: phyloseq objects and diversity metrics, one dataset at a time, with a per-dataset cache.
#
#   Step 1  scripts/0_proc_data.R <DS>         -> Out/_Rdata/cache/<DS>/ps_raw.RDS, ps_rare.RDS
#   Step 2  scripts/1_diversity_metrics.R <DS> -> Out/_Rdata/cache/<DS>/alpha_*.RDS, pcoa.RDS, pairwise_dist.RDS,
#                                                 procrustes.RDS, permanova.RDS
#   Step 3  assemble all datasets' caches      -> the usual Out/_Rdata/*.RDS files read by the figure scripts
#
# A dataset whose cache is complete isn't recomputed: adding a dataset only computes that one.
# To recompute datasets anyway (e.g. after a parsing fix, or new profiles), list them in `rerun` below.
# The diversity step is also rerun when its cache is older than the dataset's ps_rare.RDS.
#
# Adding a dataset:
#   1. HPC: run the profilers, add a line to scripts/src/datasets.tsv, then
#      bash scripts/src/merge_tables.sh <DS>   (and scripts/src/sylph_merge.sh <DS> <db> <sample tsv> per sylph database)
#   2. scripts/utils/sample_data.R: add the dataset's sample data to load_sample_data()'s list
#   3. scripts/0_config.R: add grouping_variable and dataset_names entries
#   4. add <DS> to `datasets` below and run this script (the checks below say what's missing)
#
# Run from the repo root (HPC: module load StdEnv/2023 r/4.4.0; Rscript scripts/0_run_all.R), or source() it.
library(pacman)
p_load(tidyverse, magrittr, phyloseq)

datasets <- c('P19_Saliva', 'P19_Gut', 'Moss', 'NAFLD', 'AD_Skin', 'RA_Gut', 'Bee', 'PD', 'Maple')
rerun <- c() # datasets to recompute even if their cache is complete, e.g. c('Moss')

cache_root <- 'Out/_Rdata/cache'
proc_files <- c('ps_raw.RDS', 'ps_rare.RDS')
div_files <- c('alpha_plot_data.RDS', 'alpha_wilcox.RDS', 'pcoa.RDS', 'pairwise_dist.RDS', 'procrustes.RDS', 'permanova.RDS')
dir.create('logs', showWarnings = FALSE)

source('scripts/0_config.R')
source('scripts/utils/parse_merged.R') # tool_dirs(), to check the merged tables
source('scripts/utils/sample_data.R')

# Steps 1 and 2 start a fresh R process per dataset (packages, parallel workers and memory are released after each),
# with its output in logs/run_all_<DS>_<script>.log
run_script <- function(script, ds) {
  log <- file.path('logs', paste0('run_all_', ds, '_', tools::file_path_sans_ext(basename(script)), '.log'))
  message('  running ', script, ' ', ds, ' (log: ', log, ')')
  t0 <- Sys.time()
  status <- system2(file.path(R.home('bin'), 'Rscript'), c(script, ds), stdout = log, stderr = log)
  message(sprintf('  %s %s in %.1f min', if (status == 0) 'done' else 'FAILED', script,
                  as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
  if (status != 0) message('  last lines of the log:\n    ', paste(tail(readLines(log), 15), collapse = '\n    '))
  status == 0
}


# 0. Checks before parsing a dataset that isn't cached yet ---------------------
# Returns TRUE if the dataset can be parsed; prints what's missing otherwise.
sam_data <- NULL # load_sample_data() is only called if a dataset needs parsing

check_new_dataset <- function(ds) {
  ok <- TRUE
  message('  checking ', ds, ' before parsing it:')

  # 1. scripts/src/datasets.tsv: the dataset's sample tsv (used by merge_tables.sh and sylph_merge.sh)
  tsv_lines <- readLines('scripts/src/datasets.tsv') %>% grep('^#', ., invert = TRUE, value = TRUE)
  if (!any(str_detect(tsv_lines, paste0('^', ds, '\t')))) {
    message('    MISSING: no line for ', ds, ' in scripts/src/datasets.tsv (<dataset><TAB><sample tsv>)')
    ok <- FALSE
  }

  # 2. merged tables: one per tool folder of data/<DS>/
  expected <- c()
  if (!dir.exists(file.path('data', ds))) {
    message('    MISSING: no data/', ds, '/ folder (HPC: profilers, then scripts/src/merge_tables.sh ', ds, ')')
    ok <- FALSE
  } else {
    dirs <- tool_dirs(ds)
    expected <- c(
      file.path(c(dirs$motus3, dirs$motus4, dirs$kraken, dirs$metaphlan), 'merged_profile.tsv'),
      file.path(dirs$sourmash, 'merged_tax_summary.csv.gz'),
      file.path(rep(dirs$sylph, each = 2), c('merged_species.tsv', 'merged_species_tax.tsv'))
    )
    if (length(expected) == 0) {
      message('    MISSING: no tool folder (MOTUS*, SM_*, SYLPH_*, KB*, MPA_*) in data/', ds, '/')
      ok <- FALSE
    }
  }
  missing_tables <- expected[!file.exists(file.path('data', ds, expected))]
  for (f in missing_tables) {
    if (str_detect(f, '^SYLPH_')) {
      message('    MISSING: data/', ds, '/', f, ' (HPC: bash scripts/src/sylph_merge.sh ', ds, ' ',
              str_remove(dirname(f), '^SYLPH_'), ' <sample tsv>)')
    } else {
      message('    MISSING: data/', ds, '/', f, ' (HPC: bash scripts/src/merge_tables.sh ', ds, ')')
    }
  }
  if (length(missing_tables) > 0) ok <- FALSE
  if (length(expected) > 0) message('    ', length(expected) - length(missing_tables), ' of ', length(expected), ' merged tables found')

  # 3. sample data: an entry in load_sample_data(), with sample names matching the merged tables
  if (is.null(sam_data)) sam_data <<- load_sample_data()
  if (is.null(sam_data[[ds]])) {
    message('    MISSING: no sample data for ', ds, ' in scripts/utils/sample_data.R (load_sample_data()\'s list)')
    ok <- FALSE
  } else {
    # sample names of the merged tables: the header fields of the tsv tables, query_name of the sourmash tables
    found <- expected[file.exists(file.path('data', ds, expected))]
    table_samples <- c()
    for (f in file.path('data', ds, found)) {
      if (str_ends(f, 'csv.gz')) {
        table_samples <- c(table_samples, unique(read_csv(f, col_select = 'query_name', show_col_types = FALSE)$query_name))
      } else {
        header <- readLines(f, n = 20) %>% grep('^(clade_name|mOTU|#Classification)\t', ., value = TRUE) %>% head(1)
        if (length(header) == 0) header <- readLines(f, n = 1)
        table_samples <- c(table_samples, strsplit(header, '\t')[[1]][-1])
      }
    }
    table_samples <- unique(table_samples)
    sample_ids <- rownames(sam_data[[ds]])
    n_match <- length(intersect(sample_ids, table_samples))
    message('    sample data: ', n_match, ' of ', length(sample_ids), ' samples found in the merged tables')
    if (n_match == 0 && length(found) > 0) {
      message('    MISSING: no sample of the sample data is in the merged tables: check the rownames (sample IDs) ',
              'returned for ', ds, ' by load_sample_data()')
      ok <- FALSE
    }
  }

  # 4. config: grouping variable (for the tests) and display name (for the figures). Not blocking.
  if (!ds %in% names(grouping_variable)) {
    message('    WARNING: no grouping_variable entry for ', ds, ' in scripts/0_config.R: no Wilcoxon tests or PERMANOVA')
  }
  if (!ds %in% names(dataset_names)) {
    message('    WARNING: no dataset_names entry for ', ds, ' in scripts/0_config.R: no figure label')
  }

  if (ok) message('  checks passed') else message('  ', ds, ' skipped: fix the MISSING items above')
  ok
}


# 1-2. Per dataset: parsing + rarefaction, then diversity ----------------------
for (ds in datasets) {
  message('### ', ds)
  cache_dir <- file.path(cache_root, ds)
  force <- ds %in% rerun

  # Step 1: parsing and rarefaction
  proc_cached <- all(file.exists(file.path(cache_dir, proc_files)))
  if (proc_cached && !force) {
    message('  ps_raw / ps_rare cached, skipping parsing')
  } else {
    if (force) message('  in `rerun`: parsing again')
    if (!check_new_dataset(ds)) next
    if (!run_script('scripts/0_proc_data.R', ds)) next
  }

  # Step 2: diversity metrics
  div_paths <- file.path(cache_dir, div_files)
  div_cached <- all(file.exists(div_paths))
  div_outdated <- div_cached && any(file.mtime(div_paths) < file.mtime(file.path(cache_dir, 'ps_rare.RDS')))
  if (div_cached && !div_outdated && !force) {
    message('  diversity metrics cached, skipping')
  } else {
    if (div_outdated) message('  diversity cache older than ps_rare.RDS: computing again')
    run_script('scripts/1_diversity_metrics.R', ds)
  }
}


# 3. Assemble all datasets ------------------------------------------------------
message('### Assembling')
complete <- c()
for (ds in datasets) {
  missing_files <- c(proc_files, div_files)[!file.exists(file.path(cache_root, ds, c(proc_files, div_files)))]
  if (length(missing_files) == 0) {
    complete <- c(complete, ds)
  } else {
    message('  ', ds, ' left out: missing ', paste(missing_files, collapse = ', '), ' in ', file.path(cache_root, ds))
  }
}
if (length(complete) == 0) stop('No dataset with a complete cache: nothing assembled (Out/_Rdata/ left as it was)')
message('  assembling ', length(complete), ' datasets: ', paste(complete, collapse = ', '))

read_cache <- function(ds, file) readRDS(file.path(cache_root, ds, file))

# phyloseq objects: List$Dataset$Database
ps_raw.ls <- list()
ps_rare.ls <- list()
for (ds in complete) {
  ps_raw.ls[[ds]] <- read_cache(ds, 'ps_raw.RDS')
  ps_rare.ls[[ds]] <- read_cache(ds, 'ps_rare.RDS')
}
write_rds(ps_raw.ls, 'Out/_Rdata/ps_raw.ls.RDS', compress = 'bz2')
write_rds(ps_rare.ls, 'Out/_Rdata/ps_rare.ls.RDS', compress = 'bz2')

# alpha diversity (the Wilcoxon tests in Dataset, Database, Index order, as group_split() gave them on all datasets)
alpha_div <- list()
alpha_div[['plot_data']] <- map(complete, read_cache, file = 'alpha_plot_data.RDS') %>% bind_rows()
alpha_div[['wilcox_tests']] <- map(complete, read_cache, file = 'alpha_wilcox.RDS') %>% bind_rows() %>%
  arrange(Dataset, Database, Index)
write_rds(alpha_div, 'Out/_Rdata/alpha_div.RDS', compress = 'gz')

# beta diversity
pcoa.ls <- list()
for (ds in complete) pcoa.ls[[ds]] <- read_cache(ds, 'pcoa.RDS')
write_rds(pcoa.ls, 'Out/_Rdata/pcoa_noVST.ls.RDS', compress = 'bz2')

pairwise_distances <- map(complete, read_cache, file = 'pairwise_dist.RDS') %>% bind_rows()
write_rds(pairwise_distances, 'Out/_Rdata/pairwise_dist.RDS', compress = 'bz2')

procrustes_tests <- map(complete, read_cache, file = 'procrustes.RDS') %>% bind_rows()
write_rds(procrustes_tests, 'Out/_Rdata/procrustes.RDS', compress = 'gz')

permanova.ds <- map(complete, read_cache, file = 'permanova.RDS') %>% bind_rows()
write_rds(permanova.ds, 'Out/_Rdata/permanova.ds.RDS', compress = 'gz')

# Dataset labels with n for the figures (moved from 1_diversity_metrics.R): n = samples in the Wilcoxon tests
dataset_names_n_tibble <- alpha_div[['wilcox_tests']] %>%
  mutate(n = n1 + n2) %>% # rstatix reports the two group sizes
  distinct(Dataset, n) %>%
  group_by(Dataset) %>%
  filter(n==max(n)) %>%
  left_join(tibble(
    Dataset = names(dataset_names),
    Dataset_name = unname(dataset_names)
  ), by = "Dataset") %>%
  mutate(Dataset_name_n = paste0(Dataset_name, "\n(n=",n,")" ))

# values:
dataset_names_n_lookup <- dataset_names_n_tibble$Dataset_name_n
# add standard names:
names(dataset_names_n_lookup) <- dataset_names_n_tibble$Dataset
# reorder and flush NAs
dataset_names_n_lookup %<>% .[names(dataset_names)] %>% .[!is.na(.)]
write_rds(dataset_names_n_lookup, 'Out/_Rdata/dataset_names_n_lookup_alphadiv.RDS')

message('### Written to Out/_Rdata/: ps_raw.ls, ps_rare.ls, alpha_div, pcoa_noVST.ls, pairwise_dist, procrustes, ',
        'permanova.ds, dataset_names_n_lookup_alphadiv (.RDS)')
