library(pacman)
p_load(tidyverse)

# DNA TO DNA METHODS
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
sample_subset <- class_rate %>% 
  group_by(Dataset, Sample) %>% 
  summarise(n = n()) %>% 
  filter(n > 8) # should have 10 entries per sample

# Create a label for each dataset with n = sample_count
Dataset_n_labels <- sample_subset %>% 
  select(-Sample) %>% 
  group_by(Dataset) %>% 
  summarise(n = n()) %>% 
  mutate(Dataset_n = paste0(Dataset, " (n = ", n, ")")) %>% 
  select(-n)

# Plot 
class_rate_kb %>% 
  filter(Sample %in% pull(sample_subset, Sample)
         & Dataset != 'Olive') %>% 
  ggplot(aes(x = Methodology, y = Rate, fill = Database)) +
  geom_violin(scale = 'width', 
              linewidth = 0.2, alpha = 0.8,
               position = position_dodge(width=0.6)
              ) +
  facet_grid(.~Dataset) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = c(0.25, 0.8),
    legend.background = element_rect(
      fill = "white",        # White background
      color = "black",       # Black border
      linewidth = 0.2        # Border thickness
    )) +
  scale_fill_brewer(palette = 'Set2') +
  scale_x_discrete(expand = expansion(mult = 0.2)) + # Reduce buffer between boxes and panel
  labs(x = 'Method', y = 'Read classification rate', 
       fill = 'Reference database') 

ggsave('Out/memoire/classrate_kb.pdf', bg = 'white', width = 2200, height = 1000, 
       units = 'px', dpi = 200)

 #NEXT : check by dataset type ??

class_rate_sm %>% 
  filter(Sample %in% pull(sample_subset, Sample)
         & Dataset != 'Olive') %>% 
  ggplot(aes(x = NA, y = Rate, fill = Methodology)) +
  geom_violin(scale = 'width', 
              linewidth = 0.2, alpha = 0.8,
              position = position_dodge(width=0.6)) +
  facet_grid(.~Dataset, scales = 'free', space = 'free') +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(0.1,'cm'),
    position_dodge(width = 0.7),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = c(0.25, 0.6),
    legend.background = element_rect(
      fill = "white",        # White background
      color = "black",       # Black border
      linewidth = 0.2        # Border thickness
  )) +
  scale_fill_brewer(palette = 'Pastel2')+
  scale_x_discrete(expand = expansion(mult = 0.4)) + # Reduce buffer between boxes and panel
  labs(y = 'Metagenome containment in database', 
       fill = 'Reference database') 

ggsave('Out/memoire/classrate_sm.pdf', bg = 'white', width = 2200, height = 1000, 
       units = 'px', dpi = 200)

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
) %>% select(-classified, -unclassified) %>% 
  filter(!is.na(Rate))

class_rate_markers %>% 
  filter(Sample %in% pull(sample_subset, Sample)
         & Dataset != 'Olive') %>% 
  ggplot(aes(x = NA, y = Rate, fill = Methodology)) +
  geom_violin(#scale = 'count', 
              linewidth = 0.2, alpha = 0.8) +
  facet_grid(.~Dataset, scales = 'free', space = 'free') +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(0.1,'cm'),
    position_dodge(width = 0.7),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = c(0.22, 0.6),
    legend.background = element_rect(
      fill = "white",        # White background
      color = "black",       # Black border
      linewidth = 0.2        # Border thickness
    )) +
  scale_fill_brewer(palette = 'Pastel2')+
  scale_x_discrete(expand = expansion(mult = 0.4)) + # Reduce buffer between boxes and panel
  labs(y = 'Classification rate of identified marker genes', 
       fill = 'Methodology') 

ggsave('Out/memoire/classrate_Markers.pdf', bg = 'white', width = 2200, height = 1000, 
       units = 'px', dpi = 180)
