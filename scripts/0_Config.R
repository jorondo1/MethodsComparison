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
  RA_Gut = NA
  
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
  'MOTUS' = 'mOTUs3',
  'MPA_db2019' = 'MetaPhlAn3 (2019)',
  'MPA_db2022' = 'Metaphlan4 (2022)',
  'MPA_db2023' = 'Metaphlan4 (2023)',
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

CCE_metadata <- tibble(
  Database = c("MOTUS","MPA_db2022","MPA_db2023",
               "KB10","KB45","KB90",
               "KB10_GTDB","KB45_GTDB","KB90_GTDB",
               "SM_genbank-2022.03", 'SM_RefSeq_20250528',
               "SM_gtdb-rs214-full", "SM_gtdb-rs214-rep",
               "SM_gtdb-rs220-rep","SM_gtdb-rs214-rep_MAGs"),
  MethodName = c('mOTUs', 'MPA 2022', 'MPA 2023',
                 'KB10 RefSeq', 'KB45 RefSeq', 'KB90 RefSeq',
                 'KB10 GTDB Rep.', 'KB45 GTDB Rep.', 'KB90 GTDB Rep.',
                 'SM GenBank', 'SM RefSeq',
                 'SM GTDB Full', 'SM GTDB Rep.',
                 'SM GTDB Rep.', 'SM GTDB Rep.+ MAGs'),
  MethodNameParam = c('mOTUs', 'MetaPhlAn', 'MetaPhlAn',
                      'Kraken 0.10', 'Kraken 0.45', 'Kraken 0.90',
                      'Kraken 0.10', 'Kraken 0.45', 'Kraken 0.90',
                      'Sourmash', 'Sourmash',
                      'Sourmash', 'Sourmash',
                      'Sourmash', 'Sourmash'),
  Num_species_in_db = c(25314, 30094, 36333,
                        rep(27285,3),
                        rep(113104,3),
                        62052, 27285,#Genbank has >1M genomes, but only that many bacterial and archeal species
                        85205, 85205, #GTDB full and rep have the same number of species, just more genomes
                        113104, 113211),
  plot_colour = c("#FDBF6F","green4","#3d8f58",
                  "indianred1","#c196d6","#fa817f", 
                  "indianred4", "#b41f1f", "#6A3D9A",
                  "purple3", '#87c1e0',
                  "navyblue", 'blue',
                  "#1F78B4", "skyblue3"),
  Tool = c('mOTUs3', 'MetaPhlAn4','MetaPhlAn4',
           rep('Kraken2+Bracken', 6),
           rep('Sourmash gather',6)),
  CCE_approach = c('DNA-to-Marker','DNA-to-Marker','DNA-to-Marker',
                   'DNA-to-DNA','DNA-to-DNA','DNA-to-DNA',
                   'DNA-to-DNA','DNA-to-DNA','DNA-to-DNA',
                   'DNA-to-DNA','DNA-to-DNA',
                   'DNA-to-DNA','DNA-to-DNA',
                   'DNA-to-DNA','DNA-to-DNA'),
  Taxonomy = c('Tool-specific','Tool-specific','Tool-specific',
               'NCBI','NCBI','NCBI',
               'GTDB','GTDB','GTDB',
               'NCBI','NCBI',
               'GTDB','GTDB',
               'GTDB','GTDB'),
  refdb = c('MOTUS', 'MPA', 'MPA',
            'RefSeq', 'RefSeq', 'RefSeq',
            'GTDB_220', 'GTDB_220', 'GTDB_220',
            'GenBank', 'RefSeq',
            'GTDB_214', 'GTDB_214', 
            'GTDB_220', 'GTDB_214')
)


tooldb_colours <- CCE_metadata$plot_colour
names(tooldb_colours) <- CCE_metadata$Database

tool_colours <- c(
  'mOTUs3' = "#FDBF6F", 
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
