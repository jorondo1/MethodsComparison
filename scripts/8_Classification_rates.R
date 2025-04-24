library(pacman)
p_load(tidyverse)

class_rate_kb <- read_tsv('Out/classification_rates/kraken_classification_rate.tsv', 
                       col_names = c('path','classified','unclassified','Rate')) %>% 
  separate(path, into = c('dot','Dataset', 'Tool', 'Sample', 'Report' ),
           sep = '/') %>% 
  select(-dot, -Report) %>% 
  filter(!is.na(Rate)) %>% 
  mutate(Database = case_when(str_detect(Tool, "GTDB") ~ 'GTDB rep. espèces',
                              TRUE ~ 'RefSeq'),
         Tool = str_remove(Tool, '_GTDB'))

class_rate_sm <- read_tsv('Out/classification_rates/sourmash_classification_rate.tsv',
                          col_names = c('path', 'Rate')) %>% 
  separate(path, into = c('dot', 'Dataset', 'Tool', 'Report'),
           sep = '/') %>% 
  separate(Report, into= c('Sample', NA), 
           sep = '_', extra = 'drop') %>% 
  select(-dot) %>% 
  filter(!is.na(Rate)) %>% 
  mutate(Database = case_when(str_detect(Tool, "genbank") ~ 'Genbank (NCBI)',
                              str_detect(Tool, 'rep') ~ 'GTDB rep. espèces',
                              str_detect(Tool, 'full') ~ 'GTDB complet'),
         Tool = str_remove(Tool, 'SM_')) %>% 
  filter(!is.na(Database))
  
class_rate <- rbind(
  class_rate_sm,
  class_rate_kb %>% select(-classified, -unclassified)
)

# Subset samples with both GTDB and RefSeq
sample_subset <- class_rate %>% 
  group_by(Dataset, Sample) %>% 
  summarise(n = n()) %>% 
  filter(n > 9) # should have 10 entries per sample

# Create a label for each dataset with n = sample_count
Dataset_n_labels <- sample_subset %>% 
  select(-Sample) %>% 
  group_by(Dataset) %>% 
  summarise(n = n()) %>% 
  mutate(Dataset_n = paste0(Dataset, " (n = ", n, ")")) %>% 
  select(-n)

# Plot 
class_rate_kb %>% 
  left_join(Dataset_n_labels, by = 'Dataset') %>% 
  filter(Sample %in% pull(sample_subset, Sample)) %>% 
  ggplot(aes(x = Tool, y = Rate, fill = Database)) +
  geom_boxplot(width = 0.3, linewidth = 0.2,
               position = position_dodge(width=0.5),
               outlier.size = 0.2) +
  facet_grid(.~Dataset_n) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = c(0.22, 0.8),
    legend.background = element_rect(
      fill = "white",        # White background
      color = "black",       # Black border
      linewidth = 0.2        # Border thickness
    )) +
  labs(x = 'Méthode', y = 'Taux de classification des lectures', 
       fill = 'Base de données\nde référence') 

ggsave('Out/comite2/classrate_kb.png', bg = 'white', width = 2200, height = 1000, 
       units = 'px', dpi = 200)

 #NEXT : check by dataset type ??

class_rate_sm %>% 
  left_join(Dataset_n_labels, by = 'Dataset') %>% 
  filter(Sample %in% pull(sample_subset, Sample)) %>% 
  ggplot(aes(x = NA, y = Rate, fill = Database)) +
  geom_boxplot(width = 0.3, linewidth = 0.2,
               position = position_dodge(width=0.5),
               outlier.size = 0.2) +
  facet_grid(.~Dataset_n, scales = 'free') +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = c(0.25, 0.6),
    legend.background = element_rect(
      fill = "white",        # White background
      color = "black",       # Black border
      linewidth = 0.2        # Border thickness
  )) +
  scale_fill_brewer(palette = 'Pastel2')+
  labs(x = 'Méthode', y = 'Taux de classification des lectures', 
       fill = 'Base de données\nde référence') 

ggsave('Out/comite2/classrate_sm.png', bg = 'white', width = 2200, height = 1000, 
       units = 'px', dpi = 200)
