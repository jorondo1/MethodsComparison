# Sample data (sample IDs as rownames) of each dataset, for assemble_phyloseq(). Moved from 0_proc_data.R.
# Returns list(<dataset> = data.frame). Needs tidyverse, phyloseq.
load_sample_data <- function() {
  psSalivaKB <- read_rds('data/P19_Saliva/psSalivaKB.rds')
  psFecesKB <- read_rds('data/P19_Gut/psFecesKB.rds')
  moss.ps <- readRDS(url('https://github.com/jorondo1/borealMoss/raw/main/data/R_out/mossMAGs.RDS'))

  NAFLD_meta <- read_delim('data/NAFLD/raw/ENA_report.tsv', show_col_types = FALSE) %>%
    dplyr::filter(library_strategy == 'WGS') %>%
    dplyr::select(run_accession, sample_title) %>%
    left_join(read_delim('data/NAFLD/raw/metadata.tsv', show_col_types = FALSE),
              join_by(sample_title == SampleID)) %>%
    mutate(Group = case_when(is.na(NAFLD) ~ 'Positive',
                             TRUE ~ NAFLD), .keep = 'unused') %>%
    mutate(Group = as.factor(case_when(Group == 'Negative' ~ 0,
                                       Group == 'Positive' ~ 1))) %>%
    column_to_rownames('sample_title')

  AD_skin_meta <- read_delim('data/AD_Skin/raw/ENA_report.tsv', show_col_types = FALSE) %>%
    dplyr::select(run_accession, sample_alias, run_accession) %>%
    mutate(sample_alias = str_remove(sample_alias, ' ')) %>%
    right_join(read_delim('data/AD_Skin/raw/metadata.tsv', show_col_types = FALSE),
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

  Bee_meta <- read_delim('data/Bee/raw/filereport_read_run_PRJNA685398_tsv.txt', show_col_types = FALSE) %>%
    dplyr::select(run_accession, library_name) %>%
    mutate(Group = as.factor(str_extract_all(library_name, "[A-Za-z]") %>%
                               sapply(paste, collapse = "")),
           .keep = 'unused') %>%
    column_to_rownames('run_accession')

  PD_meta <- read_delim('data/PD/raw/filereport_read_run_PRJNA834801_tsv.txt', show_col_types = FALSE) %>%
    dplyr::select(run_accession, library_name) %>%
    mutate(Group = as.factor(
      str_extract(library_name, regex(paste(c('DP', 'DC'), collapse = "|")))),
      .keep = 'unused') %>%
    filter(!is.na(Group)) %>%
    column_to_rownames('run_accession')

  # Maple: no metadata yet, only the sample IDs (data/Maple/raw/metadata.tsv: the IDs of the preprocessed reads' sample tsv)
  Maple_meta <- read_delim('data/Maple/raw/Nanodrop_omx_ME.csv', show_col_types = FALSE, delim = ";") %>%
    dplyr::select(SAMPLE_ID, Type, Provenance, traitement) %>% 
    column_to_rownames('SAMPLE_ID')

  list(
    P19_Saliva = psSalivaKB@sam_data,
    P19_Gut = psFecesKB@sam_data,
    Moss = moss.ps@sam_data,
    NAFLD = NAFLD_meta,
    AD_Skin = AD_skin_meta,
    RA_Gut = RA_meta,
    Bee = Bee_meta,
    PD = PD_meta,
    Maple = Maple_meta,
    moss.ps = moss.ps # Moss's SM_gtdb-rs214-rep_MAGs object comes from this external RDS
  )
}
