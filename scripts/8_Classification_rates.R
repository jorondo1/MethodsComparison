library(pacman)
p_load(tidyverse)

crate <- read_tsv('kraken_classification_rate.tsv', col_names = c('path','Rate')) %>% 
  separate(path, into = c('dot','Dataset', 'Tool', 'Sample', 'Report' ),
           sep = '/') %>% 
  select(-dot, -Report) %>% 
  filter(!Tool %in% c('KB20', 'KB51') & 
           !is.na(Rate)) %>% 
  mutate(Database = case_when(str_detect(Tool, "GTDB") ~ 'GTDB',
                              TRUE ~ 'RefSeq'),
         Tool = str_remove(Tool, '_GTDB'))


crate %>% 
  ggplot(aes(x = Tool, y = Rate, fill = Database)) +
  geom_boxplot() +
  facet_grid(~Dataset) +
  theme_light()

