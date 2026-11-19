if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) rstudioapi::applyTheme("Fastverse Light") # RStudio only (not Rscript on the HPC)
ggplot2::theme_set(ggplot2::theme_light())

conflicted::conflicts_prefer(
  dplyr::select, 
  dplyr::filter,
  .quiet = TRUE)

my_datasets_factorlevels <- c('P19_Saliva', 'P19_Gut', 'RA_Gut', 'AD_Skin', 'Moss', 'NAFLD')

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
  Grapevine = 'Health',   # alternative in the sample data: Compartment (paired by Plant)
  Horse = 'island_side'   # alternative in the sample data: age_class_rev
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
  Oilfield = 'Oilfield',
  Horse = 'Horse'
)

these_datasets <- names(dataset_names) # datasets to use for plots
stopifnot(!anyDuplicated(names(dataset_names)), !anyDuplicated(dataset_names))

human_datasets <- c('PD', 'RA_Gut', 'NAFLD', 'P19_Gut', 'P19_Saliva', 'AD_Skin')

# One colour per dataset (RColorBrewer Paired, pale yellow darkened)
dataset_colours <- setNames(c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C',
                              '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#E5C100', '#B15928'), these_datasets)

# Method triads compared in 2_alpha_compare.R: abundance type > reference > c(Database = label), in loop order
alpha_triads <- list(
  `Sequence abundances` = list(
    `GTDB r220` = c(KB45_GTDB = 'Kraken 0.45', SYLPH_gtdb_rs220_Seq = 'sylph', `SM_gtdb-rs220-rep` = 'sourmash'),
    RefSeq      = c(KB45 = 'Kraken 0.45', SYLPH_RefSeq_Seq = 'sylph', SM_RefSeq_20250528 = 'sourmash')),
  `Taxonomic abundances` = list(
    `Tool-specific` = c(MPA_db2025 = 'MetaPhlAn 4', MOTUS4 = 'mOTUs 4', SYLPH_gtdb_rs220_Tax = 'sylph (GTDB)')))

# Figure labels with n, e.g. "Gut PD\n(n=719)", named by dataset and in dataset_names order.
# The name comes from dataset_names above; only n (samples in the Wilcoxon tests) is read from
# dataset_names_n_lookup_alphadiv.RDS (built by 0_run_all.R), so renaming here needs no HPC rerun.
dataset_names_n <- function(sep = '\n') {
  lookup <- readRDS('Out/_Rdata/dataset_names_n_lookup_alphadiv.RDS')
  ds <- intersect(names(dataset_names), names(lookup))
  n <- sub('.*\\(n=([0-9]+)\\)$', '\\1', lookup[ds])
  setNames(paste0(dataset_names[ds], sep, '(n=', n, ')'), ds)
}

# Per-sample rates (Out/classification_rates/, from scripts/src/classification_rates.sh on the HPC), used by 6.x:
# Dataset, Database (tool dir, named like the phyloseq objects), Sample, then the file's columns
read_rates <- function(tool, columns) {
  read_tsv(paste0('Out/classification_rates/', tool, '_classification_rate.tsv'),
           col_names = c('path', columns), show_col_types = FALSE) %>%
    mutate(Dataset = str_split_i(path, '/', 3),
           Database = recode(str_split_i(path, '/', 4), MOTUS = 'MOTUS3', `SM_gtdb-rs220-reps` = 'SM_gtdb-rs220-rep'),
           Sample = str_remove(basename(path), '(_profile\\.(txt|log)|_(RefSeq|gtdb|genbank).*_gather\\.csv|_bracken_S\\.MPA\\.TXT|\\.kreport)$'),
           .before = 1, .keep = 'unused') %>%
    filter(Dataset %in% these_datasets)
}

# Post-QC samples (KneadData logs: pairs, reads = pairs x 2 + orphans), minus scripts/src/dropped_samples.tsv:
# the samples of the diversity comparisons
samples_qc <- function() {
  dropped <- read_tsv('scripts/src/dropped_samples.tsv', comment = '#', col_names = c('Dataset', 'Sample', 'reason'),
                      show_col_types = FALSE)
  read_rates('reads', c('pairs', 'reads')) %>% anti_join(dropped, by = c('Dataset', 'Sample'))
}

# Samples with a value for the grouping variable: the only ones in the Wilcoxon tests and PERMANOVA
samples_grouped <- function() {
  read_rds('Out/_Rdata/ps_raw.ls.RDS') %>%
    imap_dfr(\(l, ds) map_dfr(l, \(p) {
      sd <- as(phyloseq::sample_data(p), 'data.frame')
      tibble(Dataset = ds, Sample = rownames(sd), group = as.character(sd[[grouping_variable[[ds]]]] %||% NA))
    })) %>%
    filter(!is.na(group)) %>%
    distinct(Dataset, Sample)
}

# Dataset labels "<name> (n samples)", named by dataset
dataset_names_n_qc <- function(samples) {
  n <- count(samples, Dataset)
  setNames(paste0(dataset_names[n$Dataset], ' (', n$n, ')'), n$Dataset)
}

# Databases of the study in the rate files (sylph: one per database, not _Seq/_Tax), in display order
rate_databases <- c(MPA_db2022 = 'MetaPhlAn 2022', MPA_db2023 = 'MetaPhlAn 2023', MPA_db2025 = 'MetaPhlAn 2025',
                    MOTUS3 = 'mOTUs 3', MOTUS4 = 'mOTUs 4',
                    KB10 = 'Kraken 0.10 RefSeq', KB45 = 'Kraken 0.45 RefSeq', KB90 = 'Kraken 0.90 RefSeq',
                    KB10_GTDB = 'Kraken 0.10 GTDB', KB45_GTDB = 'Kraken 0.45 GTDB', KB90_GTDB = 'Kraken 0.90 GTDB',
                    SM_RefSeq_20250528 = 'sourmash RefSeq', `SM_gtdb-rs220-rep` = 'sourmash GTDB',
                    SYLPH_RefSeq = 'sylph RefSeq', SYLPH_gtdb_rs220 = 'sylph GTDB')

# Methods of the study: one row per profile object (Database), one column per methodological dimension.
# Replaced CCE_metadata, CCE_names, tooldb_colours and tool_colours on 2026-11-19 (legacy copy in scripts/archive/).
#   Reference: the database the tool profiles against (MetaPhlAn and mOTUs: their own marker database)
#   Param:     what distinguishes methods sharing a tool and reference (Kraken confidence, database year, version)
#   Abundance: Sequence (read counts, rarefied) or Taxonomic (proportions, not rarefied)
#   Approach:  D2M = DNA-to-marker, D2D = DNA-to-DNA
methods <- tribble(
  ~Database,              ~Tool,       ~Reference,  ~Param, ~Abundance,  ~Approach, ~colour,
  'MPA_db2022',           'MetaPhlAn', 'MetaPhlAn', '2022', 'Taxonomic', 'D2M',     'green4',
  'MPA_db2023',           'MetaPhlAn', 'MetaPhlAn', '2023', 'Taxonomic', 'D2M',     '#3d8f58',
  'MPA_db2025',           'MetaPhlAn', 'MetaPhlAn', '2025', 'Taxonomic', 'D2M',     '#00A759',
  'MOTUS3',               'mOTUs',     'mOTUs',     '3.0',    'Taxonomic', 'D2M',     'orange',
  'MOTUS4',               'mOTUs',     'mOTUs',     '4.0',    'Taxonomic', 'D2M',     '#FDBF6F',
  'KB10',                 'Kraken',    'RefSeq',    '0.10', 'Sequence',  'D2D',     'indianred1',
  'KB45',                 'Kraken',    'RefSeq',    '0.45', 'Sequence',  'D2D',     'indianred3',
  'KB90',                 'Kraken',    'RefSeq',    '0.90', 'Sequence',  'D2D',     'indianred4',
  'KB10_GTDB',            'Kraken',    'GTDB r220', '0.10', 'Sequence',  'D2D',     'purple',
  'KB45_GTDB',            'Kraken',    'GTDB r220', '0.45', 'Sequence',  'D2D',     'purple2',
  'KB90_GTDB',            'Kraken',    'GTDB r220', '0.90', 'Sequence',  'D2D',     'purple3',
  'SM_RefSeq_20250528',   'sourmash',  'RefSeq',    NA,     'Sequence',  'D2D',     '#b41f1f',
  'SM_gtdb-rs220-rep',    'sourmash',  'GTDB r220', NA,     'Sequence',  'D2D',     '#1F78B4',
  'SYLPH_RefSeq_Seq',     'sylph',     'RefSeq',    NA,     'Sequence',  'D2D',     '#8C510A',
  'SYLPH_RefSeq_Tax',     'sylph',     'RefSeq',    NA,     'Taxonomic', 'D2D',     '#BF812D',
  'SYLPH_gtdb_rs220_Seq', 'sylph',     'GTDB r220', NA,     'Sequence',  'D2D',     '#D8B365',
  'SYLPH_gtdb_rs220_Tax', 'sylph',     'GTDB r220', NA,     'Taxonomic', 'D2D',     '#F6E8C3'
) %>%
  mutate(Tool = factor(Tool, levels = c('Kraken', 'sourmash', 'sylph', 'MetaPhlAn', 'mOTUs')),
         Reference = factor(Reference, levels = c('GTDB r220', 'RefSeq', 'MetaPhlAn', 'mOTUs')),
         Abundance = factor(Abundance, levels = c('Sequence', 'Taxonomic')))
stopifnot(!anyDuplicated(methods$Database))

reference_colours <- c(`GTDB r220` = '#823D51', RefSeq = '#515A82', MetaPhlAn = '#00A759', mOTUs = '#FDBF6F')
tool_shapes <- c(Kraken = 24, sourmash = 22, sylph = 23, MetaPhlAn = 21, mOTUs = 25) # fillable shapes (21-25)

# Full label of a method, e.g. "Kraken 0.45 (GTDB r220)", "MetaPhlAn 2022", "sylph seq (RefSeq)".
# Unknown Database ids come back unchanged.
method_label <- function(db) {
  m <- methods[match(db, methods$Database), ]
  label <- case_when(
    m$Tool %in% c('MetaPhlAn', 'mOTUs') ~ paste(m$Tool, m$Param),
    m$Tool == 'sylph' ~ paste0('sylph ', if_else(m$Abundance == 'Sequence', 'seq', 'tax'), ' (', m$Reference, ')'),
    is.na(m$Param) ~ paste0(m$Tool, ' (', m$Reference, ')'),
    .default = paste0(m$Tool, ' ', m$Param, ' (', m$Reference, ')'))
  setNames(coalesce(label, db), db)
}

# Short label, for plots that already show the reference: "Kraken 0.45", "sourmash", "sylph seq", "2022" (MetaPhlAn),
# "3.0" (mOTUs)
method_short <- function(db) {
  m <- methods[match(db, methods$Database), ]
  label <- case_when(
    m$Tool %in% c('MetaPhlAn', 'mOTUs') ~ m$Param,
    m$Tool == 'sylph' ~ paste('sylph', if_else(m$Abundance == 'Sequence', 'seq', 'tax')),
    is.na(m$Param) ~ as.character(m$Tool),
    .default = paste(m$Tool, m$Param))
  setNames(coalesce(label, db), db)
}

approach_labels <- c(D2M = 'DNA-to-marker', D2D = 'DNA-to-DNA')

# Facet columns for plots with one column per method (4_alpha_tests.R, 5_beta_tests.R), for
# facet_nested(~ Approach + Reference + Method): Approach and Reference from methods, Method = method_short()
add_method_facets <- function(data) {
  data %>%
    left_join(dplyr::select(methods, Database, Approach, Reference), by = 'Database') %>%
    mutate(Approach = factor(approach_labels[Approach], levels = approach_labels),
           Method = factor(method_short(Database), levels = unique(method_short(methods$Database))))
}
# 
# tool_vars <- tibble(
#   "Aldex2" = "wi.eBH",
#   "ANCOMBC2" = "adj_p",
#   "radEmu" = "pval",
#   "corncob" = "p_fdr",
#   "DESeq2" = "padj",
#   "edgeR" = "FDR", 
#   "MaAsLin2" = "qval"
# )
# 
# DAA_metadata <- tibble(
#   DAA_tool = names(tool_vars),
#   Compositional = c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE),
#   Distribution = c('DM', 'LL', 'LL', 'BB', 'NB', 'NB', 'N'),
#   Transformation = c('CLR', 'LT', 'LT', 'NONE', 'NONE', 'NONE', 'AST'),
#   Taxon_bias = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
#   Target_ = c('RA', 'AA','AA','RA','RA','RA','RA'), # Absolute or relative abundance
#   plot_shape = c(15, 19, 17, 18, 4, 3, 6)
# )

strip_theme <- list(
  strip.text = element_text(colour = 'grey20'),
  strip.background = element_rect(
    fill = 'grey90', colour = 'grey90', linewidth = 0.1),
  panel.border = element_rect(colour = 'grey90'),
  panel.grid = element_blank(),
  panel.spacing = unit(0, "lines")   # default is 0.5, shrink gutters between panels
  
)


strip_theme_dark <- list(
  strip.text = element_text(colour = 'white'),
  strip.background = element_rect(
    fill = 'grey50', colour = 'grey50', linewidth = 0.1),
  panel.border = element_rect(colour = 'grey50'),
  panel.grid = element_blank(),
  panel.spacing = unit(0, "lines")   # default is 0.5, shrink gutters between panels
  
)
