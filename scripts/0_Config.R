rstudioapi::applyTheme("Idle Fingers") 
ggplot2::theme_set(ggplot2::theme_light())

my_datasets_factorlevels <- c('P19_Saliva', 'P19_Gut', 'RA_Gut', 'AD_Skin', 'Moss', 'NAFLD')

grouping_variable <- c(
  AD_Skin = 'Gender',
  Moss = 'Compartment',
  NAFLD = 'Group',
  P19_Gut = 'diarr',
  P19_Saliva = 'diarr',
  PD = 'Group',
  Bee = 'Group',
  Olive = NA,
  RA_Gut = 'Group'
)

dataset_names <- c(
  PD = 'Gut PD',
  RA_Gut = 'Gut RA',
  NAFLD = 'Gut NAFLD',
  P19_Gut = 'Gut COV',
  P19_Saliva = 'Saliva COV',
  AD_Skin = 'Skin AD',
  Moss = 'Moss',
  Bee = 'Bee',
  Olive = NA
)

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
  'SM_gtdb-rs214-rep_MAGs'= 'Sourmash (GTDB \n+ Novel MAGs)'
)

Hill_numbers <- c(
  'H_0' = 'Richness',
  'H_1' = 'Shannon (Hill order 1)',
  'H_2' = 'Simpson (Hill order 2)'
)

CCE_metadata <- tribble(
  ~Database,                ~plot_colour, ~MethodName,         ~MethodNameParam, ~Num_spec, ~Tool,        ~tool_family,  ~CCE_approach,   ~Taxonomy,       ~refdb,      ~short_alpha_2,      ~short_alpha_3,
  "MOTUS3",                "orange",      "mOTUs 3",           "mOTUs",              25314, "mOTUs3",     "mOTUs",       "D2M",           "Tool-specific", "MOTUS",      "mOTUs",           "3.0",
  "MOTUS4",                "#FDBF6F",     "mOTUs 4",           "mOTUs",             124296, "mOTUs4",     "mOTUs",       "D2M",           "Tool-specific", "MOTUS",      "mOTUs",           "4.0",
  "MPA_db2022",            "green4",      "MPA 2022",          "MetaPhlAn",          30094, "MetaPhlAn4", "MetaPhlAn",   "D2M",           "Tool-specific", "MPA",        "MetaPhlAn",       "2022",
  "MPA_db2023",            "#3d8f58",     "MPA 2023",          "MetaPhlAn",          36333, "MetaPhlAn4", "MetaPhlAn",   "D2M",           "Tool-specific", "MPA",        "MetaPhlAn",       "2023",
  "MPA_db2025",            "#00A759",     "MPA 2025",          "MetaPhlAn",          56335, "MetaPhlAn4", "MetaPhlAn",   "D2M",           "Tool-specific", "MPA",        "MetaPhlAn",       "2025",
  "KB10",                  "indianred1",  "KB10 RefSeq",       "Kraken 0.10",        27285, "Kraken2",    "Kraken",      "D2D",           "NCBI",          "RefSeq",     "RefSeq",  "KB 0.10",
  "KB45",                  "indianred3",  "KB45 RefSeq",       "Kraken 0.45",        27285, "Kraken2",    "Kraken",      "D2D",           "NCBI",          "RefSeq",     "RefSeq",  "KB 0.45",
  "KB90",                  "indianred4",  "KB90 RefSeq",       "Kraken 0.90",        27285, "Kraken2",    "Kraken",      "D2D",           "NCBI",          "RefSeq",     "RefSeq",  "KB 0.90",
  "KB10_GTDB",             "purple",      "KB10 GTDB",         "Kraken 0.10",       113104, "Kraken2",    "Kraken",      "D2D",       "GTDB",          "GTDB_220",   "GTDB",        "KB 0.10",
  "KB45_GTDB",             "purple2",     "KB45 GTDB",         "Kraken 0.45",       113104, "Kraken2",    "Kraken",      "D2D",       "GTDB",          "GTDB_220",   "GTDB",        "KB 0.45",
  "KB90_GTDB",             "purple3",     "KB90 GTDB",         "Kraken 0.90",       113104, "Kraken2",    "Kraken",      "D2D",       "GTDB",          "GTDB_220",   "GTDB",        "KB 0.90",
  "SM_genbank-2022.03",    "red2",        "SM GenBank",        "Sourmash",           62052, "Sourmash",   "Sourmash",    "D2D",       "NCBI",          "GenBank",    "RefSeq",  "SM",   
  "SM_RefSeq_20250528",    "#b41f1f",     "SM RefSeq",         "Sourmash",           27285, "Sourmash",   "Sourmash",    "D2D",       "NCBI",          "RefSeq",     "RefSeq",  "SM",   
  "SM_gtdb-rs214-full",    "navyblue",    "SM GTDB 214 Full",  "Sourmash",           85205, "Sourmash",   "Sourmash",    "D2D",       "GTDB",          "GTDB_214",   "GTDB",        "SM",   
  "SM_gtdb-rs214-rep",     "blue",        "SM GTDB 214",       "Sourmash",           85205, "Sourmash",   "Sourmash",    "D2D",       "GTDB",          "GTDB_214",   "GTDB",        "SM",   
  "SM_gtdb-rs220-rep",     "#1F78B4",     "SM GTDB 220",       "Sourmash",          113104, "Sourmash",   "Sourmash",    "D2D",       "GTDB",          "GTDB_220",   "GTDB",        "SM",   
  "SM_gtdb-rs214-rep_MAGs","skyblue3",    "SM GTDB 214 + MAGs","Sourmash",          113211, "Sourmash",   "Sourmash",    "D2D",       "GTDB",          "GTDB_214",   "GTDB",        "SM"
)

tooldb_colours <- CCE_metadata$plot_colour
names(tooldb_colours) <- CCE_metadata$Database

tool_colours <- c(
  'mOTUs3' = "orange", 
  'mOTUs4' = "#FDBF6F", 
  'MetaPhlAn4'= "#00A759",
  'Kraken2' = "#b41f1f",
  'Sourmash' = "#1F78B4" 
)

tool_vars <- tibble(
  "Aldex2" = "wi.eBH",
  "ANCOMBC2" = "adj_p",
  "radEmu" = "pval",
  "corncob" = "p_fdr",
  "DESeq2" = "padj",
  "edgeR" = "FDR", 
  "MaAsLin2" = "qval"
)

DAA_metadata <- tibble(
  DAA_tool = names(tool_vars),
  Compositional = c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE),
  Distribution = c('DM', 'LL', 'LL', 'BB', 'NB', 'NB', 'N'),
  Transformation = c('CLR', 'LT', 'LT', 'NONE', 'NONE', 'NONE', 'AST'),
  Taxon_bias = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  Target_ = c('RA', 'AA','AA','RA','RA','RA','RA'), # Absolute or relative abundance
  plot_shape = c(15, 19, 17, 18, 4, 3, 6)
)

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
