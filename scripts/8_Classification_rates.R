library(pacman)
p_load(tidyverse, patchwork)
theme_set(theme_light())

#########################
# DNA TO DNA METHODS ####
#########################
class_rate_kb <- read_tsv('Out/classification_rates/kraken_classification_rate.tsv', 
                          col_names = c('path','classified','unclassified','Rate')) %>% 
  separate(path, into = c(#'data',
    'dot','Dataset', 'Methodology', 'Sample', 'Report' ),
    sep = '/') %>% 
  select(-dot, -Report) %>% 
  filter(!is.na(Rate)) %>% 
  mutate(Database = case_when(str_detect(Methodology, "GTDB") ~ 'GTDB Rep.',
                              TRUE ~ 'RefSeq'),
         Methodology = str_remove(Methodology, '_GTDB'),
         Database_version = case_when(str_detect(Methodology, "GTDB") ~ '220',
                                      TRUE ~ '2024-12-28'),
         Tool = 'KB')

class_rate_sm <- read_tsv('Out/classification_rates/sourmash_classification_rate.tsv',
                          col_names = c('path', 'Rate')) %>% 
  separate(path, into = c('data','dot', 'Dataset', 'Methodology', 'Report'),
           sep = '/') %>% 
  separate(Report, into= c('Sample', NA), 
           sep = '_', extra = 'drop') %>% 
  select(-dot, -data) %>% 
  filter(!is.na(Rate)) %>% 
  mutate(Database = case_when(str_detect(Methodology, "genbank") ~ 'Genbank',
                              str_detect(Methodology, 'RefSeq') ~ 'RefSeq',
                              str_detect(Methodology, 'MAG') ~ 'GTDB Rep. + MAGs',
                              str_detect(Methodology, 'rep') ~ 'GTDB Rep.',
                              str_detect(Methodology, 'full') ~ 'GTDB Full'),
         Database_version = str_remove(Methodology, 'SM_'),
         Tool = 'SM'
  ) %>% 
  filter(!is.na(Database))

class_rate_DNA <- rbind(
  class_rate_sm,
  class_rate_kb %>% select(-classified, -unclassified)
)

# Subset samples with both GTDB and RefSeq
sample_subset <- class_rate_DNA %>% 
  group_by(Dataset, Sample) %>% 
  summarise(n = n()) %>% 
  filter(n == 11) # should have 10 entries per sample

# Create a label for each dataset with n = sample_count
Dataset_n_labels <- sample_subset %>% 
  select(-Sample) %>% 
  group_by(Dataset) %>% 
  summarise(n = n()) %>% 
  mutate(Dataset_n = paste0(Dataset, " (n = ", n, ")")) %>% 
  select(-n)

# Plot 
plot_kb <- class_rate_kb %>% 
  filter(Sample %in% pull(sample_subset, Sample)
         & Dataset != 'Olive') %>% 
  ggplot(aes(x = Methodology, y = Rate, fill = Database)) +
  geom_violin(scale = 'width', 
              linewidth = 0.2, alpha = 0.8,
              position = position_dodge(width=0.6),
  ) +
  facet_grid(.~Dataset) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  ) +
  scale_fill_brewer(palette = 'Set2') +
  scale_x_discrete(expand = expansion(mult = 0.2)) + # Reduce buffer between boxes and panel
  labs(x = 'Method', y = 'Read classification rate', 
       fill = 'Reference database',
       tag = 'A') ; plot_kb
#NEXT : check by dataset type ??

plot_sm <- class_rate_sm %>% 
  filter(Sample %in% pull(sample_subset, Sample)
         & Dataset != 'Olive') %>% 
  ggplot(aes(x = NA, y = Rate, fill = Methodology)) +
  geom_violin(scale = 'width', 
              linewidth = 0.2, alpha = 0.8,
              position = position_dodge(width=0.6)) +
  facet_grid(.~Dataset, scales = 'free', space = 'free') +
  theme(
    panel.grid.major.x = element_blank(),
    position_dodge(width = 0.7),
    axis.text.x = element_blank(),
  ) +
  scale_fill_brewer(palette = 'Pastel2')+
  scale_x_discrete(expand = expansion(mult = 0.4)) + # Reduce buffer between boxes and panel
  labs(y = 'Metagenome containment in database', 
       fill = 'Reference database',
       tag = 'B')  ; plot_sm

plot_kb / plot_sm &
  theme(    
    panel.grid.major = element_blank(),
    panel.spacing = unit(0.1,'cm'),
    # position_dodge(width = 0.7),
    strip.text = element_text(color = "black",size = 12),
    strip.background = element_rect(fill = 'grey80'),
    axis.title.x = element_blank(),
    legend.position = c(0.18, 0.6),
    legend.justification = c(0,0.5),
    legend.background = element_rect(
      fill = "white",        # White background
      color = "black",       # Black border
      linewidth = 0.2        # Border thickness
    )
  )
ggsave('Out/memoire/classrate_dnadna.pdf', bg = 'white', width = 2200, height = 2600, 
       units = 'px', dpi = 180)

############################
# DNA TO MARKER METHODS ####
############################

# Metaphlan
class_rate_mpa <- read_tsv('Out/classification_rates/mpa_classification_rate.tsv', 
                           col_names = c('path','classified','unclassified','Rate')) %>% 
  separate(path, into = c('dot','data', 'Dataset', 'Methodology', 'Sample','Report'),
           sep = '/') %>% 
  select(-dot, -Report, -data)

# MOTUS
class_rate_motus <- read_tsv('Out/classification_rates/motus_classification_rate.tsv', 
                             col_names = c('path','classified','unclassified','Rate')) %>% 
  separate(path, into = c('dot','data', 'Dataset', 'Tool', 'Report'),
           sep = '/') %>% 
  separate(Report, into = c('Sample'), sep = '_', extra = 'drop') %>% 
  select(-dot, -data, -Tool) %>% 
  mutate(Methodology = 'mOTUs 3.1.0 2023-03-28')

class_rate_markers <- rbind(
  class_rate_mpa ,
  class_rate_motus 
) %>%  filter(!is.na(Rate))

plot_marker <- class_rate_markers %>% 
  filter(Sample %in% pull(sample_subset, Sample)
         & !Dataset %in% c('Olive')
         & classified>0) %>%
  ggplot(aes(x = NA, y = Rate, fill = Methodology)) +
  geom_violin(#scale = 'count', 
    linewidth = 0.2, alpha = 0.5,
    position = position_dodge(width=0.6)) +
  facet_grid(.~Dataset, scales = 'free', space = 'free') +
  theme(
    axis.text.x = element_blank(),
  ) +
  scale_fill_brewer(palette = 'Set1', direction = -1)+
  scale_x_discrete(expand = expansion(mult = 0.4)) + # Reduce buffer between boxes and panel
  labs(y = 'Classification rate of identified marker genes', 
       fill = 'Methodology',
       tag = 'C') 

plot_kb/plot_sm/plot_marker &
  theme(    
    panel.grid.major = element_blank(),
    panel.spacing = unit(0.1,'cm'),
    # position_dodge(width = 0.7),
    strip.text = element_text(color = "black",size = 12),
    strip.background = element_rect(fill = 'grey80'),
    axis.title.x = element_blank(),
    legend.position = c(0.15, 0.6),
    legend.justification = c(0,0.5),
    legend.background = element_rect(
      fill = "white",        # White background
      color = "black",       # Black border
      linewidth = 0.2        # Border thickness
    )
  )

ggsave('Out/memoire/classrate.pdf', bg = 'white', width = 2200, height = 3000, 
       units = 'px', dpi = 180)

# Just the sourmash plot, en francais pour mémoire

db_labs_sourmash <- c(
  "SM_genbank-2022.03" = 'Genbank 2022-03', 
  'SM_RefSeq_20250528' = 'RefSeq 2024-12',
  "SM_gtdb-rs214-full" = 'GTDB 214 complet', 
  "SM_gtdb-rs214-rep" = 'GTDB 214 représentants',
  "SM_gtdb-rs220-rep" = "GTDB 220 représentants"
)
plot_palette <- RColorBrewer::brewer.pal(5, 'Pastel2')
plot_sm +
  labs(fill = 'Base de données de référence',
       y = "Proportion de k-mers du métagénome contenus dans la BD") +
  theme(legend.position = c(0.28, 0.75),
      axis.title.x = element_blank(),
      legend.background = element_rect(
        fill = "white",    
        color = "black",   
        linewidth = 0.5    
      ),
      strip.background = element_rect(fill = 'grey50'),
      strip.text.x.top = element_text(size = 12)
      )+
  scale_fill_manual(values = plot_palette,
                    labels = db_labs_sourmash)

ggsave('Out/memoire/classrate_SM_francais.pdf', 
       bg = 'white', width = 2200, height = 1400, 
       units = 'px', dpi = 180)


################################
# COMPOUND PLOT for mémoire ####
################################

combined_classrate <- rbind(
  # Kraken
  read_tsv('Out/classification_rates/kraken_classification_rate.tsv', 
           col_names = c('path','classified','unclassified','Rate')) %>% 
    separate(path, into = c(#'data',
      'dot','Dataset', 'Methodology', 'Sample', 'Report' ),
      sep = '/') %>% 
    select(Dataset, Methodology, Sample, Rate) %>% 
    filter(!is.na(Rate)),
  # Sourmash
  read_tsv('Out/classification_rates/sourmash_classification_rate.tsv',
           col_names = c('path', 'Rate')) %>% 
    separate(path, into = c('data','dot', 'Dataset', 'Methodology', 'Report'),
             sep = '/') %>% 
    separate(Report, into= c('Sample', NA), 
             sep = '_', extra = 'drop') %>% 
    select(Dataset, Methodology, Sample, Rate) %>% 
    filter(!is.na(Rate))
  )

tool_names <- c(
  "MOTUS" = "mOTUs 3.1.0",
  "MPA_db2022" = "MetaPhlAn4 vOct22",
  "MPA_db2023" = "MetaPhlAn4 vJun23",
  "KB10" = "Kraken2 0.10 RefSeq",
  "KB45" = "Kraken2 0.45 RefSeq",
  "KB90" = "Kraken2 0.90 RefSeq",
  "KB10_GTDB" = "Kraken2 0.10 GTDB 220",
  "KB45_GTDB" = "Kraken2 0.45 GTDB 220",
  "KB90_GTDB" = "Kraken2 0.90 GTDB 220",
  "SM_genbank-2022.03" = "Sourmash Genbank 2022.03", 
  'SM_RefSeq_20250528' = "Sourmash RefSeq",
  "SM_gtdb-rs214-full" = "Sourmash GTDB 214 full", 
  "SM_gtdb-rs214-rep" = "Sourmash GTDB 214",
  "SM_gtdb-rs220-rep" = "Sourmash GTDB 220"
)
tool_colours <- c(
  "MOTUS" = "#FDBF6F",
  "MPA_db2022" = "#00A759",
  "MPA_db2023" = "forestgreen",
  "KB10" = "#f28889",
  "KB45" = "#c9494b",
  "KB90" = "#a60003",
  "KB10_GTDB" = "#d77fe3",
  "KB45_GTDB" = "#a531b5",
  "KB90_GTDB" = "#560161",
  "SM_genbank-2022.03" = "#055575", 
  'SM_RefSeq_20250528' = "#3895ba",
  "SM_gtdb-rs214-full" = "#0000b8", 
  "SM_gtdb-rs214-rep" = "#383da8",
  "SM_gtdb-rs220-rep" = "#6c72eb"
)


combined_classrate %>% 
  filter(Dataset %in% c('PD', 'P19_Saliva', 'Moss') 
         & Methodology != 'SM_gtdb_rs220') %>% 
  ggplot(aes(x = Methodology, y = Rate, fill = Methodology)) +
  geom_violin(scale = 'width') +
  facet_grid(.~Dataset) +
  scale_fill_manual(values = tool_colours) +
  theme(axis.text.x = element_blank())
# mutate(
#   Tool = case_when(
#     str_detect(Methodology, 'MPA') ~ 'MetaPhlAn',
#     TRUE ~ 'mOTUs'),
#   Database = case_when(
#     str_detect(Methodology, 'MPA') ~ str_split_i(Methodology, "_", 2),
#     str_detect(Methodology, 'mOTUs') ~ str_split_i(Methodology, " ", 2))
# )


