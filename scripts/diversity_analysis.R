library(pacman)
p_load(magrittr, tidyverse, purrr, furrr, phyloseq, rstatix, foreach, doParallel)

Div_long <- read_rds('Out/Diversity_long.rds')
ps_species.ls <- read_rds("Out/ps_rare_species.ls.rds") 

add_samData <- function(df, ds, Rank, ps_list) {
  df %<>% filter(dataset == !!ds & Rank == !!Rank)
  
  # Extract sam_data from any phyloseq object of that dataset
  samData <- ps_list[[ds]][[1]]@sam_data %>% 
    as('data.frame') %>% 
    rownames_to_column('Sample')
  
  left_join(df, samData, by = 'Sample')
}

dataset <- add_samData(
  Div_long, 'Saliva','Species',ps_species.ls
  ) %>% 
  group_by(database, index) 

dataset %>%
  shapiro_test(value) %>%
  filter(p<0.05)

# Check consistency of testing diversity difference between male/female
test_results <- dataset %>% 
  wilcox_test(value ~ sex) %>% # conservative
  add_significance() %>% 
  select(database, index, p.signif) %>%
  mutate(p.signif = factor(p.signif, levels = c('ns','*','**','***', '****'))) # reorder factors

# Plot every database x index combination
test_results %>% # add the p-values to dataset
  left_join(dataset, join_by('database', 'index')) %>% 
  ggplot(aes(x = sex, y = value, colour = p.signif)) +
  scale_color_discrete() +
  geom_boxplot(linewidth = 0.3) + 
  facet_grid(index~database, scales = 'free_y') +
  theme(
    panel.background = element_rect(fill = "grey95")  # Sets a lighter grey background
  )

ggsave('Out/alpha_saliva_sex.pdf', bg = 'white', 
       width = 1800, height = 2400, units = 'px', dpi = 240)
