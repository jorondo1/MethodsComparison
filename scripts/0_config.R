# Shared configuration, sourced by every analysis script.
#
# Contents
#   1. Setup
#   2. Datasets
#   3. Methods
#      3.1 Methods table
#      3.2 Labels
#   4. Aesthetics
#      4.1 Colours and shapes
#      4.2 Themes
#   5. Figure settings
#   6. Samples
#   7. Classification rates

# 1. Setup ----------------------------------------------------------------

# RStudio only (not Rscript on the HPC)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  rstudioapi::applyTheme("Fastverse Light")
}
ggplot2::theme_set(ggplot2::theme_light())

conflicted::conflicts_prefer(
  dplyr::select,
  dplyr::filter,
  .quiet = TRUE)

# 2. Datasets -------------------------------------------------------------

# Variable tested in each dataset (Wilcoxon, PERMANOVA)
grouping_variable <- c(
  AD_Skin = 'Gender',
  Moss = 'Compartment',
  NAFLD = 'Group',
  P19_Gut = 'diarr',
  P19_Saliva = 'diarr',
  PD = 'Group',
  Bee = 'Group',
  RA_Gut = 'Group',
  Maple = 'Provenance',
  Oilfield = 'Oilfield',
  Grapevine = 'Health',   # alternative: Compartment (paired by Plant)
  Horse = 'island_side'   # alternative: age_class_rev
)

# Display names; this order is the dataset order in every figure
dataset_names <- c(
  PD = 'Gut PD',
  RA_Gut = 'Gut RA',
  NAFLD = 'Gut NAFLD',
  P19_Gut = 'Gut COV',
  P19_Saliva = 'Saliva COV',
  AD_Skin = 'Skin AD',
  Bee = 'Bee',
  Moss = 'Moss',
  Maple = 'Maple',
  Grapevine = 'Grapevine',
  Oilfield = 'Oil field',
  Horse = 'Horse'
)

these_datasets <- names(dataset_names)

human_datasets <- c('PD', 'RA_Gut', 'NAFLD', 'P19_Gut', 'P19_Saliva', 'AD_Skin')

# 3. Methods --------------------------------------------------------------

## 3.1 Methods table ----
# One row per profile object (Database)
#   Reference: database profiled against (MetaPhlAn, mOTUs: their own markers)
#   Param:     what tells apart methods sharing tool and reference
#   Abundance: Sequence (read counts, rarefied) or Taxonomic (proportions)
#   Approach:  D2M = DNA-to-marker, D2D = DNA-to-DNA
#   Short:     abbreviation for crowded axes (RS = RefSeq)
#   colour:    per-method colour (10_benchmark_plot.R)
methods <- tribble(
  ~Database,              ~Short,        ~Tool,            ~Reference,  ~Param, ~Abundance,  ~Approach, ~colour,
  'MPA_db2022',           'MPA 2022',    'MetaPhlAn',      'MetaPhlAn', '2022', 'Taxonomic', 'D2M',     'green4',
  'MPA_db2023',           'MPA 2023',    'MetaPhlAn',      'MetaPhlAn', '2023', 'Taxonomic', 'D2M',     '#3d8f58',
  'MPA_db2025',           'MPA 2025',    'MetaPhlAn',      'MetaPhlAn', '2025', 'Taxonomic', 'D2M',     '#00A759',
  'MOTUS3',               'mOTUs 3.0',   'mOTUs',          'mOTUs',     '3.0',  'Taxonomic', 'D2M',     'orange',
  'MOTUS4',               'mOTUs 4.0',   'mOTUs',          'mOTUs',     '4.0',  'Taxonomic', 'D2M',     '#FDBF6F',
  'KB10',                 'KB10 RS',     'Kraken-Bracken', 'RefSeq',    '0.10', 'Sequence',  'D2D',     'indianred1',
  'KB45',                 'KB45 RS',     'Kraken-Bracken', 'RefSeq',    '0.45', 'Sequence',  'D2D',     'indianred3',
  'KB90',                 'KB90 RS',     'Kraken-Bracken', 'RefSeq',    '0.90', 'Sequence',  'D2D',     'indianred4',
  'KB10_GTDB',            'KB10 GTDB',   'Kraken-Bracken', 'GTDB r220', '0.10', 'Sequence',  'D2D',     'purple',
  'KB45_GTDB',            'KB45 GTDB',   'Kraken-Bracken', 'GTDB r220', '0.45', 'Sequence',  'D2D',     'purple2',
  'KB90_GTDB',            'KB90 GTDB',   'Kraken-Bracken', 'GTDB r220', '0.90', 'Sequence',  'D2D',     'purple3',
  'SM_RefSeq_20250528',   'SM RS',       'sourmash',       'RefSeq',    NA,     'Sequence',  'D2D',     '#b41f1f',
  'SM_gtdb-rs220-rep',    'SM GTDB',     'sourmash',       'GTDB r220', NA,     'Sequence',  'D2D',     '#1F78B4',
  'SYLPH_RefSeq_Seq',     'SY Seq RS',   'sylph',          'RefSeq',    NA,     'Sequence',  'D2D',     '#8C510A',
  'SYLPH_RefSeq_Tax',     'SY Tax RS',   'sylph',          'RefSeq',    NA,     'Taxonomic', 'D2D',     '#BF812D',
  'SYLPH_gtdb_rs220_Seq', 'SY Seq GTDB', 'sylph',          'GTDB r220', NA,     'Sequence',  'D2D',     '#D8B365',
  'SYLPH_gtdb_rs220_Tax', 'SY Tax GTDB', 'sylph',          'GTDB r220', NA,     'Taxonomic', 'D2D',     '#F6E8C3'
) %>%
  mutate(
    Tool = factor(
      Tool,
      levels = c('Kraken-Bracken', 'sourmash', 'sylph', 'MetaPhlAn', 'mOTUs')),
    Reference = factor(
      Reference,
      levels = c('GTDB r220', 'RefSeq', 'MetaPhlAn', 'mOTUs')),
    Abundance = factor(
      Abundance,
      levels = c('Sequence', 'Taxonomic')),
    # GTDB, RefSeq, or the tool's own marker database
    RefClass = factor(
      if_else(Reference %in% c('MetaPhlAn', 'mOTUs'), 'Tool-specific', as.character(Reference)),
      levels = c('GTDB r220', 'RefSeq', 'Tool-specific')))

stopifnot(
  !anyDuplicated(methods$Database),
  !anyDuplicated(methods$Short))

## 3.2 Labels ----
# e.g. labels = method_abbrev in a discrete scale
method_abbrev <- setNames(methods$Short, methods$Database)

# "Kraken-Bracken 0.45 (GTDB r220)", "MetaPhlAn 2022", "sylph seq (RefSeq)";
# unknown ids come back unchanged
method_label <- function(db) {
  m <- methods[match(db, methods$Database), ]
  label <- case_when(
    m$Tool %in% c('MetaPhlAn', 'mOTUs') ~ paste(m$Tool, m$Param),
    m$Tool == 'sylph' ~ paste0(
      'sylph ', if_else(m$Abundance == 'Sequence', 'seq', 'tax'),
      ' (', m$Reference, ')'),
    is.na(m$Param) ~ paste0(m$Tool, ' (', m$Reference, ')'),
    .default = paste0(m$Tool, ' ', m$Param, ' (', m$Reference, ')'))
  setNames(coalesce(label, db), db)
}

approach_labels <- c(
  D2M = 'DNA-to-marker',
  D2D = 'DNA-to-DNA')

# 4. Aesthetics -----------------------------------------------------------

## 4.1 Colours and shapes ----
# Method points (Figs 3-5): colour = tool, shape = reference class
# Okabe-Ito: colour-blind safe, still visible as light tints
tool_colours <- c(
  `Kraken-Bracken` = '#0072B2',
  sourmash = '#D55E00',
  sylph = '#CC79A7',
  MetaPhlAn = '#009E73',
  mOTUs = '#E69F00')

# fillable shapes (21-25)
refclass_shapes <- c(
  `GTDB r220` = 24,
  RefSeq = 22,
  `Tool-specific` = 21)

# Strips of 22_compare_beta.R
refclass_colours <- c(
  `GTDB r220` = '#823D51',
  RefSeq = '#515A82',
  `Tool-specific` = '#45733A')

# Neutral: hues are taken by tool_colours
approach_colours <- c(
  `DNA-to-marker` = '#C8C8C8',
  `DNA-to-DNA` = '#4D4D4D')

dataset_colours <- c(
  PD = '#E16A86',
  RA_Gut = '#00AD9A',
  NAFLD = '#B88A00',
  P19_Gut = '#E063B6',
  P19_Saliva = '#50A315',
  AD_Skin = '#9183E6',
  # non-human
  Bee = '#D37A4A',
  Moss = '#00A8C3',
  Maple = '#C86DD7',
  Grapevine = '#909800',
  Oilfield = '#00AB66',
  Horse = '#009ADE')

## 4.2 Themes ----
strip_theme <- list(
  strip.text = element_text(colour = 'grey20'),
  strip.background = element_rect(
    fill = 'grey90',
    colour = 'grey90',
    linewidth = 0.1),
  panel.border = element_rect(colour = 'grey90'),
  panel.grid = element_blank(),
  panel.spacing = unit(0, "lines"))

# 5. Figure settings ------------------------------------------------------

# Alpha indices shown, with display names
these_indices <- c(
  Richness = 'Richness',
  Shannon = 'Shannon',
  Hill_2 = 'Inverse Simpson',
  Tail = 'Tail')

# One colour per index (Okabe-Ito; 1.5.alpha_rank_consistency)
index_colours <- c(
  Richness = '#000000',
  Shannon = '#E69F00',
  `Inverse Simpson` = '#56B4E9',
  Tail = '#009E73')

# Triads of 21_compare_alpha.R: abundance type > reference > c(Database = label);
# the loop repeats the first method (sylph, common to all triads) as last axis
alpha_triads <- list(
  `Sequence abundances` = list(
    `GTDB r220` = c(
      SYLPH_gtdb_rs220_Seq = 'sylph',
      KB45_GTDB = 'Kraken-Bracken 0.45',
      `SM_gtdb-rs220-rep` = 'sourmash'),
    RefSeq = c(
      SYLPH_RefSeq_Seq = 'sylph',
      KB45 = 'Kraken-Bracken 0.45',
      SM_RefSeq_20250528 = 'sourmash')),
  `Taxonomic abundances` = list(
    `Tool-specific` = c(
      SYLPH_gtdb_rs220_Tax = 'sylph (GTDB)',
      MPA_db2025 = 'MetaPhlAn 4',
      MOTUS4 = 'mOTUs 4')))

# 6. Samples --------------------------------------------------------------

# "Gut PD\n(n=719)", in dataset_names order; n = samples in the Wilcoxon tests
# (dataset_names_n_lookup_alphadiv.RDS, from 0_run_all.R)
dataset_names_n <- function(sep = '\n') {
  lookup <- readRDS('Out/_Rdata/dataset_names_n_lookup_alphadiv.RDS')
  ds <- intersect(names(dataset_names), names(lookup))
  n <- sub('.*\\(n=([0-9]+)\\)$', '\\1', lookup[ds])
  setNames(paste0(dataset_names[ds], sep, '(n=', n, ')'), ds)
}

# Post-QC samples (KneadData logs), minus scripts/src/dropped_samples.tsv
samples_qc <- function() {
  dropped <- read_tsv(
    'scripts/src/dropped_samples.tsv',
    comment = '#',
    col_names = c('Dataset', 'Sample', 'reason'),
    show_col_types = FALSE)
  read_rates('reads', c('pairs', 'reads')) %>%
    anti_join(dropped, by = c('Dataset', 'Sample'))
}

# Samples with a grouping value: those in the Wilcoxon tests and PERMANOVA
samples_grouped <- function() {
  read_rds('Out/_Rdata/ps_raw.ls.RDS') %>%
    imap_dfr(\(l, ds) map_dfr(l, \(p) {
      sd <- as(phyloseq::sample_data(p), 'data.frame')
      tibble(
        Dataset = ds,
        Sample = rownames(sd),
        group = as.character(sd[[grouping_variable[[ds]]]] %||% NA))
    })) %>%
    filter(!is.na(group)) %>%
    distinct(Dataset, Sample)
}

# "<name> (n)", named by dataset
dataset_names_n_qc <- function(samples) {
  n <- count(samples, Dataset)
  setNames(paste0(dataset_names[n$Dataset], ' (', n$n, ')'), n$Dataset)
}

# 7. Classification rates -------------------------------------------------

# Out/classification_rates/<tool>_classification_rate.tsv (HPC):
# returns Dataset, Database, Sample, then the file's columns
read_rates <- function(tool, columns) {
  read_tsv(
    paste0('Out/classification_rates/', tool, '_classification_rate.tsv'),
    col_names = c('path', columns),
    show_col_types = FALSE) %>%
    mutate(
      Dataset = str_split_i(path, '/', 3),
      # tool dir, renamed like the phyloseq objects
      Database = recode(
        str_split_i(path, '/', 4),
        MOTUS = 'MOTUS3',
        `SM_gtdb-rs220-reps` = 'SM_gtdb-rs220-rep'),
      Sample = str_remove(
        basename(path),
        '(_profile\\.(txt|log)|_(RefSeq|gtdb|genbank).*_gather\\.csv|_bracken_S\\.MPA\\.TXT|\\.kreport)$'),
      .before = 1,
      .keep = 'unused') %>%
    filter(Dataset %in% these_datasets)
}

# Databases in the rate files (sylph: one per database), in display order
rate_databases <- c(
  MPA_db2022 = 'MetaPhlAn 2022',
  MPA_db2023 = 'MetaPhlAn 2023',
  MPA_db2025 = 'MetaPhlAn 2025',
  MOTUS3 = 'mOTUs 3',
  MOTUS4 = 'mOTUs 4',
  KB10 = 'Kraken-Bracken 0.10 RefSeq',
  KB45 = 'Kraken-Bracken 0.45 RefSeq',
  KB90 = 'Kraken-Bracken 0.90 RefSeq',
  KB10_GTDB = 'Kraken-Bracken 0.10 GTDB',
  KB45_GTDB = 'Kraken-Bracken 0.45 GTDB',
  KB90_GTDB = 'Kraken-Bracken 0.90 GTDB',
  SM_RefSeq_20250528 = 'sourmash RefSeq',
  `SM_gtdb-rs220-rep` = 'sourmash GTDB',
  SYLPH_RefSeq = 'sylph RefSeq',
  SYLPH_gtdb_rs220 = 'sylph GTDB')
