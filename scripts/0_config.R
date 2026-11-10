if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) rstudioapi::applyTheme("Idle Fingers") # RStudio only (not Rscript on the HPC)
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

CCE_names <- c(
  'MOTUS3' = 'mOTUs3',
  'MOTUS4' = 'mOTUs4',
  'MPA_db2022' = 'Metaphlan4 (2022)',
  'MPA_db2023' = 'Metaphlan4 (2023)',
  'MPA_db2025' = 'MetaPhlAn4 (2025)',
  'KB10' = 'Kraken 0.10 (RefSeq)',
  'KB45' = 'Kraken 0.45 (RefSeq)',
  'KB90' = 'Kraken 0.90 (RefSeq)',
  'KB10_GTDB' = 'Kraken 0.10 (GTDB 220)',
  'KB45_GTDB' = 'Kraken 0.45 (GTDB 220)',
  'KB90_GTDB' = 'Kraken 0.90 (GTDB 220)',
  'SM_genbank-2022.03' = 'Sourmash (Genbank)',
  'SM_RefSeq_20250528' = 'Sourmash (RefSeq)',
  'SM_gtdb-rs214-full' = 'Sourmash (GTDB 214 full)',
  'SM_gtdb-rs214-rep' = 'Sourmash (GTDB 214)',
  'SM_gtdb-rs220-rep'= 'Sourmash (GTDB 220)',
  'SM_gtdb-rs214-rep_MAGs'= 'Sourmash (GTDB \n+ Novel MAGs)',
  'SYLPH_RefSeq_Seq' = 'sylph (RefSeq, sequence)',
  'SYLPH_RefSeq_Tax' = 'sylph (RefSeq, taxonomic)',
  'SYLPH_gtdb_rs220_Seq' = 'sylph (GTDB 220, sequence)',
  'SYLPH_gtdb_rs220_Tax' = 'sylph (GTDB 220, taxonomic)'
)

Hill_numbers <- c(
  'H_0' = 'Richness',
  'H_1' = 'Shannon (Hill order 1)',
  'H_2' = 'Simpson (Hill order 2)'
)

CCE_metadata <- tribble(
  ~Database,                ~plot_colour, ~MethodName,         ~MethodNameParam, ~Num_spec, ~Tool,        ~tool_family,  ~CCE_approach,   ~Taxonomy,       ~refdb,      ~short_alpha_2,      ~short_alpha_3,
  "MOTUS3",                "orange",      "mOTUs 3",           "mOTUs",              25314, "mOTUs3",     "mOTUs",       "D2M",       "Tool-specific", "MOTUS",      "mOTUs",           "3.0",
  "MOTUS4",                "#FDBF6F",     "mOTUs 4",           "mOTUs",             124296, "mOTUs4",     "mOTUs",       "D2M",       "Tool-specific", "MOTUS",      "mOTUs",           "4.0",
  "MPA_db2022",            "green4",      "MPA 2022",          "MetaPhlAn",          30094, "MetaPhlAn4", "MetaPhlAn",   "D2M",       "Tool-specific", "MPA",        "MetaPhlAn",       "2022",
  "MPA_db2023",            "#3d8f58",     "MPA 2023",          "MetaPhlAn",          36333, "MetaPhlAn4", "MetaPhlAn",   "D2M",       "Tool-specific", "MPA",        "MetaPhlAn",       "2023",
  "MPA_db2025",            "#00A759",     "MPA 2025",          "MetaPhlAn",          56335, "MetaPhlAn4", "MetaPhlAn",   "D2M",       "Tool-specific", "MPA",        "MetaPhlAn",       "2025",
  "KB10",                  "indianred1",  "KB10 RefSeq",       "Kraken 0.10",        27285, "Kraken2",    "Kraken",      "D2D",       "NCBI",          "RefSeq",     "RefSeq",  "KB 0.10",
  "KB45",                  "indianred3",  "KB45 RefSeq",       "Kraken 0.45",        27285, "Kraken2",    "Kraken",      "D2D",       "NCBI",          "RefSeq",     "RefSeq",  "KB 0.45",
  "KB90",                  "indianred4",  "KB90 RefSeq",       "Kraken 0.90",        27285, "Kraken2",    "Kraken",      "D2D",       "NCBI",          "RefSeq",     "RefSeq",  "KB 0.90",
  "KB10_GTDB",             "purple",      "KB10 GTDB",         "Kraken 0.10",       113104, "Kraken2",    "Kraken",      "D2D",       "GTDB",          "GTDB_220",   "GTDB",        "KB 0.10",
  "KB45_GTDB",             "purple2",     "KB45 GTDB",         "Kraken 0.45",       113104, "Kraken2",    "Kraken",      "D2D",       "GTDB",          "GTDB_220",   "GTDB",        "KB 0.45",
  "KB90_GTDB",             "purple3",     "KB90 GTDB",         "Kraken 0.90",       113104, "Kraken2",    "Kraken",      "D2D",       "GTDB",          "GTDB_220",   "GTDB",        "KB 0.90",
  "SM_genbank-2022.03",    "red2",        "SM GenBank",        "Sourmash",           62052, "Sourmash",   "Sourmash",    "D2D",       "NCBI",          "GenBank",    "RefSeq",  "SM",   
  "SM_RefSeq_20250528",    "#b41f1f",     "SM RefSeq",         "Sourmash",           27285, "Sourmash",   "Sourmash",    "D2D",       "NCBI",          "RefSeq",     "RefSeq",  "SM",   
  "SM_gtdb-rs214-full",    "navyblue",    "SM GTDB 214 Full",  "Sourmash",           85205, "Sourmash",   "Sourmash",    "D2D",       "GTDB",          "GTDB_214",   "GTDB",        "SM",   
  "SM_gtdb-rs214-rep",     "blue",        "SM GTDB 214",       "Sourmash",           85205, "Sourmash",   "Sourmash",    "D2D",       "GTDB",          "GTDB_214",   "GTDB",        "SM",   
  "SM_gtdb-rs220-rep",     "#1F78B4",     "SM GTDB 220",       "Sourmash",          113104, "Sourmash",   "Sourmash",    "D2D",       "GTDB",          "GTDB_220",   "GTDB",        "SM",   
  "SM_gtdb-rs214-rep_MAGs","skyblue3",    "SM GTDB 214 + MAGs","Sourmash",          113211, "Sourmash",   "Sourmash",    "D2D",       "GTDB",          "GTDB_214",   "GTDB",        "SM",
  # Num_spec: same genome sets as SM_RefSeq_20250528 (RefSeq 20250528) and the GTDB 220 entries above (r220 representatives)
  # _Seq: estimated read counts (rarefied); _Tax: taxonomic abundance, proportions (not rarefied, like MetaPhlAn)
  "SYLPH_RefSeq_Seq",      "#8C510A",     "sylph RefSeq seq",  "sylph",              27285, "sylph",      "sylph",       "D2D",       "NCBI",          "RefSeq",     "RefSeq",  "sylph seq",
  "SYLPH_RefSeq_Tax",      "#BF812D",     "sylph RefSeq tax",  "sylph",              27285, "sylph",      "sylph",       "D2D",       "NCBI",          "RefSeq",     "RefSeq",  "sylph tax",
  "SYLPH_gtdb_rs220_Seq",  "#D8B365",     "sylph GTDB 220 seq","sylph",             113104, "sylph",      "sylph",       "D2D",       "GTDB",          "GTDB_220",   "GTDB",        "sylph seq",
  "SYLPH_gtdb_rs220_Tax",  "#F6E8C3",     "sylph GTDB 220 tax","sylph",             113104, "sylph",      "sylph",       "D2D",       "GTDB",          "GTDB_220",   "GTDB",        "sylph tax"
)

tooldb_colours <- CCE_metadata$plot_colour
names(tooldb_colours) <- CCE_metadata$Database

tool_colours <- c(
  'mOTUs3' = "orange", 
  'mOTUs4' = "#FDBF6F", 
  'MetaPhlAn4'= "#00A759",
  'Kraken2' = "#b41f1f",
  'Sourmash' = "#1F78B4",
  'sylph' = "#8C510A"
)
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
