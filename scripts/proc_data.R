library(pacman)
p_load(phyloseq, tidyverse, magrittr, doParallel)

# functions
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/community_functions.R'))
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/rarefy_even_depth2.R'))
source('scripts/myFunctions.R')

# Import External data
objectsToImport <- c("psSalivaKB", "psFecesKB")
for (i in objectsToImport) {assign(i,readRDS(
  paste0("/Users/jorondo/Library/CloudStorage/OneDrive-USherbrooke/Projets/PROVID19/objects/",i,".rds")))}
moss.ps <-readRDS(url('https://github.com/jorondo1/borealMoss/raw/main/data/R_out/mossMAGs.RDS'))

# Function to call parsing functions for each dataset,
# results in a list containing 1 phyloseq object by dataset by tool (database)
# with general structure List$Dataset$Database.ps
meta_parsing <- function(dsName, samData, filtering = FALSE) {
  ps <- list()
  # Metaphlan  
  for (db in c('MPA_db2022', 'MPA_db2023')) {
    ps[[db]] <- parse_MPA(
      MPA_files = paste0(dsName,'/', db, '/*/*profile.txt'),
      column_names = c('Taxonomy', 'NCBI','Abundance', 'Void')) %>% 
      assemble_phyloseq(samData, filtering = filtering)
  }
  
  # Kraken-bracken (using default headers from parse_MPA function)
  ps[['Bracken51']] <- parse_MPA(
    MPA_files = paste0(dsName,'/Bracken51/*/*_bracken/*_bracken_S.MPA.TXT')) %>% 
    assemble_phyloseq(samData, filtering = filtering)
  
    # MOTUS
  ps[['MOTUS']] <- parse_MPA(
    MPA_files = paste0(dsName,"/MOTUS/*_profile.txt"), 
    column_names = c('mOTU', 'Taxonomy', 'NCBI', 'Abundance')) %>% 
    assemble_phyloseq(samData, filtering = filtering)
  
  # Sourmash
  ps[['SM_genbank_202203']] <- left_join(
    parse_SM(paste0(dsName,'/SM_genbank_202203/*_genbank-2022.03_gather.csv')),
    parse_genbank_lineages(paste0(dsName,'/SM_genbank_202203/genbank-2022.03_lineages.csv')),
    by = 'genome'
    ) %>% species_glom() %>%
    assemble_phyloseq(samData, filtering = filtering)
  
  for (db in c(#'SM_gtdb_rs214_full',
               'SM_gtdb_rs214_rep')) {
    ps[[db]] <- left_join(
      parse_SM(paste0(dsName,'/', db, '/*_gather.csv')),
      parse_GTDB_lineages(paste0(dsName,'/', db, '/', db, '_lineages.csv')),
      by = 'genome'
    ) %>% assemble_phyloseq(samData, filtering = filtering)
  }
  return(ps)
}

ps.ls <- list()
ps.ls[['Saliva']] <- meta_parsing('Saliva', psSalivaKB@sam_data)
ps.ls[['Feces']] <- meta_parsing('Feces', psFecesKB@sam_data)
ps.ls[['Moss']] <- meta_parsing('Moss', moss.ps@sam_data)
ps.ls[['Moss']][['SM_gtdb_rs214_rep_MAGs']] <- moss.ps 
ps.ls$Moss$MPA_db2022 <- NULL
ps.ls$Moss$MPA_db2023 <- NULL
ps.ls$Moss$MOTUS <- NULL

# Prevalence+Abundance filtering, currently hardcoded in filter_low_prevalence()
ps_filt.ls <- list()
ps_filt.ls[['Saliva']] <- meta_parsing('Saliva', psSalivaKB@sam_data, filtering = TRUE)
ps_filt.ls[['Feces']] <- meta_parsing('Feces', psFecesKB@sam_data, filtering = TRUE)
ps_filt.ls[['Moss']] <- meta_parsing('Moss', moss.ps@sam_data, filtering = TRUE)
ps_filt.ls[['Moss']][['SM_gtdb_rs214_rep_MAGs']] <- moss.ps%>% filter_low_prevalence()
ps_filt.ls$Moss$MPA_db2022 <- NULL
ps_filt.ls$Moss$MPA_db2023 <- NULL
ps_filt.ls$Moss$MOTUS <- NULL


# Rarefy all datasets
ps_filt_rare.ls <- lapply(ps_filt.ls, function(sublist) {
  lapply(sublist, rarefy_even_depth2, rngseed = 1234)
})

ps_rare.ls <- lapply(ps.ls, function(sublist) {
  lapply(sublist, rarefy_even_depth2, rngseed = 1234)
})

write_rds(ps_filt.ls, "Out/ps_filt.ls.rds")
write_rds(ps.ls, "Out/ps.ls.rds")
write_rds(ps_filt_rare.ls, "Out/ps_filt_rare.ls.rds")
write_rds(ps_rare.ls, "Out/ps_rare.ls.rds")
