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
  AD_Skin = 'AD Skin',
  Moss = 'Moss',
  NAFLD = 'NAFLD Gut',
  P19_Gut = 'P19 Gut',
  P19_Saliva = 'P19 Saliva',
  PD = 'PD Gut',
  Bee = 'Bee',
  Olive = NA,
  RA_Gut = NA
)

CCE_names <- c(
  'MOTUS3' = 'mOTUs3',
  'MOTUS4' = 'mOTUs4',
  'MPA_db2022' = 'Metaphlan4 (2022)',
  'MPA_db2023' = 'Metaphlan4 (2023)',
  'MPA_db2025' = 'MetaPhlAn3 (2025)',
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
  ~Database,                ~plot_colour, ~MethodName,         ~MethodNameParam, ~Num_spec, ~Tool,                ~tool_family,  ~CCE_approach,      ~Taxonomy,       ~refdb,      ~short_alpha_2,      ~short_alpha_3,
  "MOTUS3",                "orange",      "mOTUs 3",           "mOTUs",              25314, "mOTUs3",             "mOTUs",       "DNA-to-Marker",    "Tool-specific", "MOTUS",      "mOTUs",           "3.0",
  "MOTUS4",                "#FDBF6F",     "mOTUs 4",           "mOTUs",             124296, "mOTUs4",             "mOTUs",       "DNA-to-Marker",    "Tool-specific", "MOTUS",      "mOTUs",           "4.0",
  "MPA_db2022",            "green4",      "MPA 2022",          "MetaPhlAn",          30094, "MetaPhlAn4",         "MetaPhlAn",   "DNA-to-Marker",    "Tool-specific", "MPA",        "MetaPhlAn",       "2022-12",
  "MPA_db2023",            "#3d8f58",     "MPA 2023",          "MetaPhlAn",          36333, "MetaPhlAn4",         "MetaPhlAn",   "DNA-to-Marker",    "Tool-specific", "MPA",        "MetaPhlAn",       "2023-07",
  "MPA_db2025",            "palegreen3",  "MPA 2025",          "MetaPhlAn",          56335, "MetaPhlAn4",         "MetaPhlAn",   "DNA-to-Marker",    "Tool-specific", "MPA",        "MetaPhlAn",       "2025-03",
  "KB10",                  "indianred1",  "KB10 RefSeq",       "Kraken 0.10",        27285, "Kraken2+Bracken",    "Kraken",      "DNA-to-DNA",       "NCBI",          "RefSeq",     "RefSeq 2024-12",  "Kraken 0.10",
  "KB45",                  "#c196d6",     "KB45 RefSeq",       "Kraken 0.45",        27285, "Kraken2+Bracken",    "Kraken",      "DNA-to-DNA",       "NCBI",          "RefSeq",     "RefSeq 2024-12",  "Kraken 0.45",
  "KB90",                  "#fa817f",     "KB90 RefSeq",       "Kraken 0.90",        27285, "Kraken2+Bracken",    "Kraken",      "DNA-to-DNA",       "NCBI",          "RefSeq",     "RefSeq 2024-12",  "Kraken 0.90",
  "KB10_GTDB",             "indianred4",  "KB10 GTDB Rep.",    "Kraken 0.10",       113104, "Kraken2+Bracken",    "Kraken",      "DNA-to-DNA",       "GTDB",          "GTDB_220",   "GTDB 220",        "Kraken 0.10",
  "KB45_GTDB",             "#b41f1f",     "KB45 GTDB Rep.",    "Kraken 0.45",       113104, "Kraken2+Bracken",    "Kraken",      "DNA-to-DNA",       "GTDB",          "GTDB_220",   "GTDB 220",        "Kraken 0.45",
  "KB90_GTDB",             "#6A3D9A",     "KB90 GTDB Rep.",    "Kraken 0.90",       113104, "Kraken2+Bracken",    "Kraken",      "DNA-to-DNA",       "GTDB",          "GTDB_220",   "GTDB 220",        "Kraken 0.90",
  "SM_genbank-2022.03",    "purple3",     "SM GenBank",        "Sourmash",           62052, "Sourmash gather",    "Sourmash",    "DNA-to-DNA",       "NCBI",          "GenBank",    "RefSeq 2024-12",  "Sourmash",   
  "SM_RefSeq_20250528",    "#87c1e0",     "SM RefSeq",         "Sourmash",           27285, "Sourmash gather",    "Sourmash",    "DNA-to-DNA",       "NCBI",          "RefSeq",     "RefSeq 2024-12",  "Sourmash",   
  "SM_gtdb-rs214-full",    "navyblue",    "SM GTDB Full",      "Sourmash",           85205, "Sourmash gather",    "Sourmash",    "DNA-to-DNA",       "GTDB",          "GTDB_214",   "GTDB 214",        "Sourmash",   
  "SM_gtdb-rs214-rep",     "blue",        "SM GTDB Rep.",      "Sourmash",           85205, "Sourmash gather",    "Sourmash",    "DNA-to-DNA",       "GTDB",          "GTDB_214",   "GTDB 214",        "Sourmash",   
  "SM_gtdb-rs220-rep",     "#1F78B4",     "SM GTDB Rep.",      "Sourmash",          113104, "Sourmash gather",    "Sourmash",    "DNA-to-DNA",       "GTDB",          "GTDB_220",   "GTDB 220",        "Sourmash",   
  "SM_gtdb-rs214-rep_MAGs","skyblue3",    "SM GTDB Rep.+ MAGs","Sourmash",          113211, "Sourmash gather",    "Sourmash",    "DNA-to-DNA",       "GTDB",          "GTDB_214",   "GTDB 214",        "Sourmash"
)

tooldb_colours <- CCE_metadata$plot_colour
names(tooldb_colours) <- CCE_metadata$Database

tool_colours <- c(
  'mOTUs3' = "orange", 
  'mOTUs4' = "#FDBF6F", 
  'MetaPhlAn4'= "#00A759",
  'Kraken2+Bracken' = "#b41f1f",
  'Sourmash gather' = "#1F78B4" 
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
