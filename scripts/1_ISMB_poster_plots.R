library(pacman)
p_load(
  # syntax
  magrittr, tidyverse, purrr, furrr, foreach, doParallel, parallel,
  # metrics and stats
  phyloseq, DESeq2, vegan, rstatix,
  # plotting :
  ggbeeswarm2, patchwork, grid, ggh4x)

ps.ls <- read_rds('Out/ps_full.ls.rds')
ps_rare.ls <- read_rds('Out/ps_rare.ls.rds')





Div_long <- read_rds('Out/Diversity_long.rds')

## Maybe one comparative diversity plot using lines only?
# The most variable lines could be thicker...?
Div_long %>% 
  filter(#dataset %in% c('Moss', 'P19_Saliva') & 
           database != 'KB05' & 
           index == 'H_1' &
           Rank == 'Species') %>% 
  # filter(Rank != 'Family' & index == 'H_0' & dataset != 'Feces') %>% 
  ggplot(aes(x = database, y = value, fill = database)) +
  geom_violin(size = 0.3) + 
  geom_line(aes(group = Sample), alpha=0.3, linewidth = 0.1) + 
  facet_nested(cols = vars(dataset), rows=vars(index), scales = 'free') +
  theme_light() + 
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values = tool_colours)
