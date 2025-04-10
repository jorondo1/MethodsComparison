library(pacman)
p_load(tidyverse)

class_rate <- read_tsv('kraken_classification_rate.tsv', col_names = c('path','Rate')) %>% 
  separate(path, into = c('dot','Dataset', 'Tool', 'Sample', 'Report' ),
           sep = '/') %>% 
  select(-dot, -Report) %>% 
  filter(!Tool %in% c('KB20', 'KB51') & 
           !is.na(Rate)) %>% 
  mutate(Database = case_when(str_detect(Tool, "GTDB") ~ 'GTDB',
                              TRUE ~ 'RefSeq'),
         Tool = str_remove(Tool, '_GTDB'))

# Subset samples with both GTDB and RefSeq
sample_subset <- class_rate %>% 
  group_by(Dataset, Sample) %>% 
  summarise(n = n()) %>% 
  filter(n == 6) # should have 6 entries per sample

# Create a label for each dataset with n = sample_count
Dataset_n_labels <- sample_subset %>% 
  select(-Sample) %>% 
  group_by(Dataset) %>% 
  summarise(n = n()) %>% 
  mutate(Dataset_n = paste0(Dataset, " (n = ", n, ")")) %>% 
  select(-n)

# Plot 
class_rate %>% 
  left_join(Dataset_n_labels, by = 'Dataset') %>% 
  filter(Sample %in% pull(sample_subset, Sample)) %>% 
  ggplot(aes(x = Tool, y = Rate, fill = Database)) +
  geom_boxplot(width = 0.5, #alpha = 0.7,
               position = position_dodge(width=0.5)) +
  facet_grid(~Dataset_n) +
  theme_light()

