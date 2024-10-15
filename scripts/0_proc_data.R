library(pacman)
p_load(phyloseq, tidyverse, magrittr, doParallel)

# functions
source('scripts/myFunctions.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/tax_glom2.R')

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
  
  # Metaphlan  
  for (db in c('MPA_db2019','MPA_db2022', 'MPA_db2023')) {
    message(paste('Parsing', db, '...'))
    ps[[db]] <- parse_MPA(
      MPA_files = paste0(dsName,'/', db, '/*/*_profile.txt'),
      column_names = c('Taxonomy', 'NCBI','Abundance', 'Void'),
      convert_to_counts = TRUE) %>% 
      assemble_phyloseq(samData)
  }
  
  # Sourmash
  message(paste('Parsing', 'SM_genbank-2022.03', '...'))
  ps[['SM_genbank-2022.03']] <- left_join(
    parse_SM(paste0(dsName,'/SM_genbank-2022.03/*_genbank-2022.03_gather.csv')),
    parse_genbank_lineages(paste0(dsName,'/SM_genbank-2022.03/SM_genbank-2022.03_lineages.csv')),
    by = 'genome'
  ) %>% species_glom() %>%
    assemble_phyloseq(samData)
  
  for (db in c('SM_gtdb-rs214-full',
               'SM_gtdb-rs214-rep')) {
    message(paste('Parsing', db, '...'))
    
    ps[[db]] <- left_join(
      parse_SM(paste0(dsName,'/', db, '/*_gather.csv')),
      parse_GTDB_lineages(paste0(dsName,'/', db, '/', db, '_lineages.csv')),
      by = 'genome'
    ) %>% species_glom() %>%
      assemble_phyloseq(samData)
  }

  # Kraken-bracken (using default headers from parse_MPA function)
  for (db in c('KB20', 'KB51')) {
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

# Retrieve NAFLD metadata
NAFLD_meta <- read_delim('NAFLD/raw/ENA_report.tsv') %>% 
  dplyr::filter(library_strategy == 'WGS') %>% 
  dplyr::select(run_accession, sample_title) %>% 
  left_join(read_delim('NAFLD/raw/metadata.tsv'), 
            join_by(sample_title == SampleID)) %>% 
  mutate(Group = case_when(is.na(NAFLD) ~ 'Positive',
                           TRUE ~ NAFLD), .keep = 'unused') %>% 
  mutate(Group = case_when(Group == 'Negative' ~ 0,
                           Group == 'Positive' ~ 1)) %>% 
  column_to_rownames('sample_title')

AD_skin_meta <- read_delim('AD_Skin/raw/ENA_report.tsv') %>% 
  dplyr::select(run_accession, sample_alias, run_accession) %>% 
  mutate(sample_alias = str_remove(sample_alias, ' ')) %>% 
  right_join(read_delim('AD_Skin/raw/metadata.tsv'),
            join_by(sample_alias == SampleID)) %>% 
  dplyr::select(-`Sample Location`) %>% 
  mutate(BGA = `Birth Gestational Age (weeks)`, 
         Group = case_when(Group == 'control' ~ 0,
                           Group == 'eczema' ~ 1),
         .keep = 'unused') %>% 
  column_to_rownames('run_accession')

# Full phyloseq objects 
ps_raw.ls <- list()
ps_raw.ls[['P19_Saliva']] <- meta_parsing('P19_Saliva', psSalivaKB@sam_data)
ps_raw.ls[['P19_Gut']] <- meta_parsing('P19_Gut', psFecesKB@sam_data)
ps_raw.ls[['Moss']] <- meta_parsing('Moss', moss.ps@sam_data)
ps_raw.ls[['NAFLD']] <- meta_parsing('NAFLD', NAFLD_meta)
ps_raw.ls[['AD_Skin']] <- meta_parsing('AD_Skin', AD_skin_meta)
ps_raw.ls[['Moss']][['SM_gtdb-rs214-rep_MAGs']] <- moss.ps 
ps_raw.ls$Moss$MPA_db2022 <- NULL
ps_raw.ls$Moss$MPA_db2023 <- NULL
ps_raw.ls$Moss$MOTUS <- NULL

# Prevalence+Abundance filtering, currently hardcoded in filter_low_prevalence()
ps_filt.ls <- lapply(ps_raw.ls, function(ds) {
  lapply(ds, filter_low_prevalence, minPrev = 0.10, minAbund = 0)
})

ps_full.ls <- list()
ps_full.ls[['Species']] <- ps_filt.ls
ps_full.ls[['Genus']] <- lapply(ps_filt.ls, function(ds) {
  lapply(ds, function(db) {
    tax_glom2(db, taxrank = "Genus") 
  })
})

ps_full.ls[['Family']] <- lapply(ps_filt.ls, function(ds) {
  lapply(ds, function(db) {
    tax_glom2(db, taxrank = "Family") 
  })
})

ps_rare.ls <- list()
ps_rare.ls[['Species']] <- lapply(ps_raw.ls, function(ds) {
  lapply(ds, function(db) {
    rarefy_even_depth2(db, rngseed = 1234)
  })
})

ps_rare.ls[['Genus']] <- lapply(ps_raw.ls, function(ds) {
  lapply(ds, function(db) {
    tax_glom2(db, taxrank = "Family") %>% 
      rarefy_even_depth2(rngseed = 1234)
  })
})

ps_rare.ls[['Family']] <- lapply(ps_raw.ls, function(ds) {
  lapply(ds, function(db) {
    tax_glom2(db, taxrank = "Family") %>% 
      rarefy_even_depth2(rngseed = 1234)
  })
})


write_rds(ps_full.ls, "Out/ps_full.ls.rds")
#write_rds(ps_filt.ls, "Out/ps_filt.ls.rds")
write_rds(ps_rare.ls, "Out/ps_rare.ls.rds")


# ps_rare.ls <- lapply(ps_raw.ls, function(sublist) {
#   lapply(sublist, rarefy_even_depth2, rngseed = 1234)
# })
# 
# # Genus level phyloseq objects, with prevalence filtering
# ps_rare_genus.ls <- lapply(ps_raw.ls, function(ds) {
#   lapply(ds, function(db) {
#     tax_glom2(db, taxrank = "Genus") %>% 
#       rarefy_even_depth2(rngseed = 1234, verbose = TRUE)
#     })
# })
# 
# ps_rare_family.ls <- lapply(ps_raw.ls, function(ds) {
#   lapply(ds, function(db) {
#     tax_glom2(db, taxrank = "Family") %>% 
#       rarefy_even_depth2(rngseed = 1234)
#   })
# })
