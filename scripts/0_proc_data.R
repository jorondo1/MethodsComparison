library(pacman)
p_load(phyloseq,mgx.tools, # unloadNamespace("mgx.tools"); devtools::install_github("jorondo1/mgx.tools", force = TRUE); library(mgx.tools)
       tidyverse, magrittr, doParallel, furrr)

dir.create('Out/_Rdata', recursive = TRUE, showWarnings = FALSE)

# Readers for the merged tables (scripts/src/merge_tables.sh on the HPC) and meta_parsing():
# one phyloseq object per database, List$Dataset$Database.ps
source('scripts/utils/parse_merged.R')
# Sample data of each dataset (+ the external Moss MAGs object)
source('scripts/utils/sample_data.R')
sam_data <- load_sample_data()

# Full phyloseq objects 
ps_raw.ls <- list()
ps_raw.ls[['P19_Saliva']] <- meta_parsing('P19_Saliva', sam_data$P19_Saliva)
ps_raw.ls[['P19_Gut']] <- meta_parsing('P19_Gut', sam_data$P19_Gut)
ps_raw.ls[['Moss']] <- meta_parsing('Moss', sam_data$Moss)
ps_raw.ls[['Moss']][['SM_gtdb-rs214-rep_MAGs']] <- sam_data$moss.ps
ps_raw.ls$Moss$MPA_db2022 <- NULL
ps_raw.ls$Moss$MPA_db2023 <- NULL
ps_raw.ls$Moss$KB90 <- NULL
ps_raw.ls$Moss$MOTUS <- NULL
ps_raw.ls[['NAFLD']] <- meta_parsing('NAFLD', sam_data$NAFLD)
ps_raw.ls[['AD_Skin']] <- meta_parsing('AD_Skin', sam_data$AD_Skin)
ps_raw.ls[['RA_Gut']] <- meta_parsing('RA_Gut', sam_data$RA_Gut)
ps_raw.ls[['Bee']] <- meta_parsing('Bee', sam_data$Bee)
ps_raw.ls[['PD']] <- meta_parsing('PD', sam_data$PD)

write_rds(ps_raw.ls, "Out/_Rdata/ps_raw.ls.RDS", compress = 'bz2')

ps_rare.ls <- imap(ps_raw.ls, function(ds, dataset_name) {
  imap(ds, function(db, database_name) {
    if (str_detect(database_name, "^MPA_")) {
      message(paste("Skipping rarefaction for:", dataset_name, database_name))
      return(db)
    } else {
      message(paste("Rarefying:",dataset_name,  database_name))
      mgx.tools::rarefy_even_depth2(
        db, rngseed = 1234, 
        verbose = TRUE, ncores = 7) %>% 
        return()
    }
  })
})

write_rds(ps_rare.ls, "Out/_Rdata/ps_rare.ls.RDS", compress = 'bz2')



# For DA project
# # Prevalence+Abundance filtering, currently hardcoded in filter_low_prevalence()
# ps_filt.ls <- lapply(ps_raw.ls, function(ds) {
#   lapply(ds, filter_low_prevalence, minPrev = 0.10, minAbund = 0)
# })
# write_rds(ps_filt.ls, "Out/_Rdata/ps_filt.ls.RDS", compress = 'bz2')
