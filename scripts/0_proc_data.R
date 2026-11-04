# Phyloseq objects of ONE dataset, from its merged tables -> Out/_Rdata/cache/<DS>/ps_raw.RDS and ps_rare.RDS
# Usually run by scripts/0_run_all.R, which then assembles all datasets into Out/_Rdata/ps_raw.ls.RDS / ps_rare.ls.RDS.
# By hand: Rscript scripts/0_proc_data.R <dataset>   (or set DS <- '<dataset>' and source() this file)
pacman::p_load(phyloseq,mgx.tools, # unloadNamespace("mgx.tools"); devtools::install_github("jorondo1/mgx.tools", force = TRUE); library(mgx.tools)
       tidyverse, magrittr, doParallel, furrr)

if (length(commandArgs(trailingOnly = TRUE)) > 0) DS <- commandArgs(trailingOnly = TRUE)[1]
if (!exists('DS')) stop('No dataset given: Rscript scripts/0_proc_data.R <dataset>')
cache_dir <- file.path('Out/_Rdata/cache', DS)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
message('### ', DS, ': parsing and rarefaction -> ', cache_dir)

# Cores for the rarefaction: the Slurm allocation on the HPC, else all but one
N_CORES <- as.integer(Sys.getenv('SLURM_CPUS_ON_NODE', parallel::detectCores() - 1))

# Readers for the merged tables (scripts/src/merge_tables.sh on the HPC), meta_parsing() and dataset_fixes()
source('scripts/utils/parse_merged.R')
# Sample data of each dataset (+ the external Moss MAGs object)
source('scripts/utils/sample_data.R')
sam_data <- load_sample_data()

# Full phyloseq objects: one per database, ps_raw$<database>
ps_raw <- meta_parsing(DS, sam_data[[DS]])
ps_raw <- dataset_fixes(DS, ps_raw, sam_data)

write_rds(ps_raw, file.path(cache_dir, 'ps_raw.RDS'), compress = 'gz')

ps_rare <- imap(ps_raw, function(db, database_name) {
  if (str_detect(database_name, "^MPA_|^SYLPH_.*_Tax$")) { # proportions
    message(paste("Skipping rarefaction for:", DS, database_name))
    return(db)
  } else {
    message(paste("Rarefying:", DS, database_name))
    # doRNG: the same result whatever the number of cores
    mgx.tools::rarefy_even_depth2(
      db, rngseed = 1234,
      verbose = TRUE, ncores = N_CORES) %>%
      return()
  }
})

write_rds(ps_rare, file.path(cache_dir, 'ps_rare.RDS'), compress = 'gz')
message('### ', DS, ': done')



# For DA project
# # Prevalence+Abundance filtering, currently hardcoded in filter_low_prevalence()
# ps_filt.ls <- lapply(ps_raw.ls, function(ds) {
#   lapply(ds, filter_low_prevalence, minPrev = 0.10, minAbund = 0)
# })
# write_rds(ps_filt.ls, "Out/_Rdata/ps_filt.ls.RDS", compress = 'bz2')
