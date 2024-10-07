library(pacman)
p_load(phyloseq, tidyverse, magrittr, doParallel)

# functions
source('scripts/myFunctions.R')

# Import External data
objectsToImport <- c("psSalivaKB", "psFecesKB")
for (i in objectsToImport) {assign(i,readRDS(
  paste0("/Users/jorondo/Library/CloudStorage/OneDrive-USherbrooke/Projets/PROVID19/objects/",i,".rds")))}
moss.ps <-readRDS(url('https://github.com/jorondo1/borealMoss/raw/main/data/R_out/mossMAGs.RDS'))

# Function to call parsing functions for each dataset,
# results in a list containing 1 phyloseq object by dataset by tool (database)
# with general structure List$Dataset$Database.ps
meta_parsing <- function(dsName, samData) {
  ps <- list()
  # Sourmash
  message(paste('Parsing', 'SM_genbank-2022.03', '...'))
  ps[['SM_genbank-2022.03']] <- left_join(
    parse_SM(paste0(dsName,'/SM_genbank-2022.03/*_genbank-2022.03_gather.csv')),
    parse_genbank_lineages(paste0(dsName,'/SM_genbank-2022.03/SM_genbank-2022.03_lineages.csv')),
    by = 'genome'
  ) %>% species_glom() %>%
    assemble_phyloseq(samData)
  
  for (db in c('SM_gtdb_rs214_full',
               'SM_gtdb_rs214_rep')) {
    message(paste('Parsing', db, '...'))
    
    ps[[db]] <- left_join(
      parse_SM(paste0(dsName,'/', db, '/*_gather.csv')),
      parse_GTDB_lineages(paste0(dsName,'/', db, '/', db, '_lineages.csv')),
      by = 'genome'
    ) %>% species_glom() %>%
      assemble_phyloseq(samData)
  }

  # Metaphlan  
  for (db in c('MPA_db2022', 'MPA_db2023')) {
    message(paste('Parsing', db, '...'))
    ps[[db]] <- parse_MPA(
      MPA_files = paste0(dsName,'/', db, '/*/*profile.txt'),
      column_names = c('Taxonomy', 'NCBI','Abundance', 'Void')) %>% 
      assemble_phyloseq(samData)
  }
  
  # Kraken-bracken (using default headers from parse_MPA function)
  for (db in c('KB05','KB20', 'KB51')) {
    message(paste('Parsing', db, '...'))
    ps[[db]] <- parse_MPA(
    MPA_files = paste0(dsName,'/', db, '/*/*_bracken/*_bracken_S.MPA.TXT')) %>% 
    assemble_phyloseq(samData)
  }
    # MOTUS
  message(paste('Parsing', 'MOTUS', '...'))
  ps[['MOTUS']] <- parse_MPA(
    MPA_files = paste0(dsName,"/MOTUS/*_profile.txt"), 
    column_names = c('mOTU', 'Taxonomy', 'NCBI', 'Abundance')) %>% 
    assemble_phyloseq(samData)
  
  return(ps)
}

# Retrieve RA_Gut metadata
RA_samdata <- read.delim('RA_Gut/raw/metadata.tsv') %>% 
  select(sra_ftp, sample_alias) %>% 
  mutate(Sample = str_extract(basename(sra_ftp), "SRR\\d+"),
         Group = str_extract(sample_alias, "[A-Za-z]+"),   # Extracts the letters part
         Num_ID = str_extract(sample_alias, "\\d+"),
         .keep = 'unused') %>% 
  column_to_rownames('Sample')

# Full phyloseq objects 
ps_raw.ls <- list()
ps_raw.ls[['P19_Saliva']] <- meta_parsing('P19_Saliva', psSalivaKB@sam_data)
ps_raw.ls[['P19_Gut']] <- meta_parsing('P19_Gut', psFecesKB@sam_data)
ps_raw.ls[['Moss']] <- meta_parsing('Moss', moss.ps@sam_data)
ps_raw.ls[['RA_Gut']] <- meta_parsing('RA_Gut', RA_samdata)
ps_raw.ls[['Moss']][['SM_gtdb_rs214_rep_MAGs']] <- moss.ps 
ps_raw.ls$Moss$MPA_db2022 <- NULL
ps_raw.ls$Moss$MPA_db2023 <- NULL
ps_raw.ls$Moss$MOTUS <- NULL

# Prevalence+Abundance filtering, currently hardcoded in filter_low_prevalence()
ps_filt.ls <- lapply(ps_raw.ls, function(ds) {
  lapply(ds, filter_low_prevalence)
})

ps_rare.ls <- lapply(ps_filt.ls, function(sublist) {
  lapply(sublist, rarefy_even_depth2, rngseed = 1234)
})

# Genus level phyloseq objects, with prevalence filtering
ps_rare_genus.ls <- lapply(ps_filt.ls, function(ds) {
  lapply(ds, function(db) {
    tax_glom(db, taxrank = "Genus") %>% 
      rarefy_even_depth2(rngseed = 1234)
    })
})

ps_rare_family.ls <- lapply(ps_filt.ls, function(ds) {
  lapply(ds, function(db) {
    tax_glom(db, taxrank = "Family") %>% 
      rarefy_even_depth2(rngseed = 1234)
  })
})

write_rds(ps_raw.ls, "Out/ps_raw.ls.rds")
write_rds(ps_filt.ls, "Out/ps_filt.ls.rds")
write_rds(ps_rare.ls, "Out/ps_rare_species.ls.rds")
write_rds(ps_rare_genus.ls, "Out/ps_rare_genus.ls.rds")
write_rds(ps_rare_family.ls, "Out/ps_rare_family.ls.rds")

