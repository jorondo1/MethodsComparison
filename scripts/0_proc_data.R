library(pacman)
p_load(phyloseq, tidyverse, magrittr, doParallel, furrr)

# functions
source('scripts/myFunctions.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/tax_glom2.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/rarefy_even_depth2.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R')

# Import External data
psSalivaKB <- read_rds('data/P19_Saliva/psSalivaKB.rds')
psFecesKB <- read_rds('data/P19_Gut/psFecesKB.rds')
moss.ps <-readRDS(url('https://github.com/jorondo1/borealMoss/raw/main/data/R_out/mossMAGs.RDS'))

# Function to call parsing functions for each dataset,
# results in a list containing 1 phyloseq object by dataset by tool (database)
# with general structure List$Dataset$Database.ps
meta_parsing <- function(dsName, samData) {
  ps <- list()
  
  # SOURMASH #####################
  SM_dirs <- list.dirs(file.path('data',dsName), recursive = FALSE) %>% 
    .[grep("/SM_[^/]*$",.)] %>% basename
  
  for (db in SM_dirs) {
    message(paste('Parsing', db, '...'))
    lineage_path <- file.path('data', dsName, db, paste0(db, '_lineages.csv'))
    
    ps[[db]] <- left_join(
      parse_SM(file.path('data', dsName, db, '*_gather.csv')),
      if(str_detect(db, 'gtdb')) {
        parse_GTDB_lineages(lineage_path)
      } else {
        parse_genbank_lineages(lineage_path)
      },
      by = 'genome'
    ) %>% species_glom() %>%
      assemble_phyloseq(samData)
  }
  
  # MOTUS
  message(paste('Parsing', 'MOTUS', '...'))
  ps[['MOTUS']] <- parse_MPA(
    MPA_files = file.path('data', dsName,"MOTUS/*_profile.txt"), 
    column_names = c('mOTU', 'Taxonomy', 'NCBI', 'Abundance'),
    mOTUs_data = TRUE) %>% 
    assemble_phyloseq(samData)
  
  # Kraken-bracken (using default headers from parse_MPA function)
  kbdirs <- list.dirs(file.path('data',dsName), recursive = FALSE) %>% 
    .[grep("/KB[^/]*$", .)] %>% basename # List all KB dirs
  
  for (db in kbdirs) {
    message(paste('Parsing', db, '...'))
    ps[[db]] <- parse_MPA(
      MPA_files = file.path('data', dsName, db, '*/*_bracken/*_bracken_S.MPA.TXT')) %>% 
      assemble_phyloseq(samData)
  }
  

  # METAPHLAN ######################
  mpadirs <- list.dirs(file.path('data',dsName), recursive = FALSE) %>% 
    .[grep("/MPA_[^/]*$", .)] %>% basename
  
  for (db in mpadirs) {
    message(paste('Parsing', db, '...'))
    ps[[db]] <- parse_MPA(
      MPA_files = file.path('data', dsName, db, '*/*_profile.txt'),
      column_names = c('Taxonomy', 'NCBI','Abundance', 'Void'),
      convert_to_counts = TRUE) %>% 
      assemble_phyloseq(samData)
  }
  return(ps)
}

# Retrieve NAFLD metadata
NAFLD_meta <- read_delim('data/NAFLD/raw/ENA_report.tsv') %>% 
  dplyr::filter(library_strategy == 'WGS') %>% 
  dplyr::select(run_accession, sample_title) %>% 
  left_join(read_delim('data/NAFLD/raw/metadata.tsv'), 
            join_by(sample_title == SampleID)) %>% 
  mutate(Group = case_when(is.na(NAFLD) ~ 'Positive',
                           TRUE ~ NAFLD), .keep = 'unused') %>% 
  mutate(Group = as.factor(case_when(Group == 'Negative' ~ 0,
                                     Group == 'Positive' ~ 1))) %>% 
  column_to_rownames('sample_title')

AD_skin_meta <- read_delim('data/AD_Skin/raw/ENA_report.tsv') %>% 
  dplyr::select(run_accession, sample_alias, run_accession) %>% 
  mutate(sample_alias = str_remove(sample_alias, ' ')) %>% 
  right_join(read_delim('data/AD_Skin/raw/metadata.tsv'),
             join_by(sample_alias == SampleID)) %>% 
  dplyr::select(-`Sample Location`) %>% 
  mutate(BGA = `Birth Gestational Age (weeks)`, 
         Group = case_when(Group == 'control' ~ 0,
                           Group == 'eczema' ~ 1),
         .keep = 'unused') %>% 
  column_to_rownames('run_accession')

RA_meta <- read.delim('data/RA_Gut/raw/metadata.tsv') %>% 
  tibble %>% 
  transmute(Sample = run_accession, 
            Group = as.factor(str_extract(sample_alias, "^[A-Za-z]+"))) %>% 
  column_to_rownames('Sample')

Bee_meta <- read_delim('data/Bee/raw/filereport_read_run_PRJNA685398_tsv.txt') %>% 
  select(run_accession, library_name) %>% 
  mutate(Group = as.factor(str_extract_all(library_name, "[A-Za-z]") %>% 
                             sapply(paste, collapse = "")),
         .keep = 'unused') %>% 
  column_to_rownames('run_accession')

Olive_meta <- read_delim('data/Olive/raw/filereport_read_run_PRJNA629675_tsv.txt') %>% 
  select(run_accession, library_name) %>% 
  mutate(Group = as.factor(
    str_extract(library_name, regex(paste(c('Kal', 'FS'), collapse = "|")))),
    .keep = 'unused') %>% 
  column_to_rownames('run_accession')

PD_meta <- read_delim('data/PD/raw/filereport_read_run_PRJNA834801_tsv.txt') %>% 
  select(run_accession, library_name) %>% 
  mutate(Group = as.factor(
    str_extract(library_name, regex(paste(c('DP', 'DC'), collapse = "|")))),
    .keep = 'unused') %>% 
  filter(!is.na(Group)) %>% 
  column_to_rownames('run_accession')

# Full phyloseq objects 
ps_raw.ls <- list()
ps_raw.ls[['P19_Saliva']] <- meta_parsing('P19_Saliva', psSalivaKB@sam_data)
ps_raw.ls[['P19_Gut']] <- meta_parsing('P19_Gut', psFecesKB@sam_data)
ps_raw.ls[['Moss']] <- meta_parsing('Moss', moss.ps@sam_data)
ps_raw.ls[['Moss']][['SM_gtdb-rs214-rep_MAGs']] <- moss.ps 
ps_raw.ls$Moss$MPA_db2022 <- NULL
ps_raw.ls$Moss$MPA_db2023 <- NULL
ps_raw.ls$Moss$KB90 <- NULL
ps_raw.ls$Moss$MOTUS <- NULL
ps_raw.ls[['NAFLD']] <- meta_parsing('NAFLD', NAFLD_meta)
ps_raw.ls[['AD_Skin']] <- meta_parsing('AD_Skin', AD_skin_meta)
#ps_raw.ls[['RA_Gut']] <- meta_parsing('RA_Gut', RA_meta)
ps_raw.ls[['Bee']] <- meta_parsing('Bee', Bee_meta)
ps_raw.ls[['Olive']] <- meta_parsing('Olive', Olive_meta)
ps_raw.ls[['PD']] <- meta_parsing('PD', PD_meta)

write_rds(ps_raw.ls, "Out/ps_raw.ls.RDS")

# Prevalence+Abundance filtering, currently hardcoded in filter_low_prevalence()

ps_filt.ls <- lapply(ps_raw.ls, function(ds) {
  lapply(ds, filter_low_prevalence, minPrev = 0.10, minAbund = 0)
})
write_rds(ps_filt.ls, "Out/ps_filt.ls.RDS")

ps_rare.ls <- list()
ps_rare.ls <- lapply(ps_raw.ls, function(ds) {
  lapply(ds, function(db) {
    rarefy_even_depth2(db, rngseed = 1234, 
                       verbose = TRUE, ncores = 6)
  })
})

write_rds(ps_rare.ls, "Out/ps_rare.ls.RDS")


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
