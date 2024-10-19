library(pacman)
p_load(
  # syntax
  magrittr, tidyverse, purrr, furrr, foreach, doParallel, parallel,
  # metrics and stats
  phyloseq, DESeq2, vegan, rstatix,
  # plotting :
  ggbeeswarm2, patchwork, grid, ggh4x, ggthemes)

ps.ls <- read_rds('Out/ps_full.ls.rds')
ps_rare.ls <- read_rds('Out/ps_rare.ls.rds')
Div_long <- read_rds('Out/Diversity_long_full.rds') %>% 
  mutate(index_div = case_when(index == 'H_1' ~ log(value),
                               index == 'H_2' ~ 1-1/value, 
                               TRUE ~ value))

## Maybe one comparative diversity plot using lines only?
# Find most variable lines per group
set.seed(37); most_variable <- Div_long %>% 
  select(index, Rank, dataset, Sample) %>%
  unique %>% 
  mutate(random = runif(n(), min = 0, max = 1)) %>% 
  filter(random > 0.94) %>% 
    mutate(alpha = 1, linewidth = 0.3) 
  
# The most variable lines could be thicker...?
Div_long %>% 
  filter(dataset %in% c('Moss') & 
           index == 'H_1' &
           Rank == 'Species') %>% 
  left_join(most_variable, join_by('index', 'Rank', 'dataset', 'Sample')) %>% 
  mutate(alpha = coalesce(alpha, 0.2),
         linewidth = coalesce(linewidth, 0.3)) %>% 
  ggplot(aes(x = database, y = index_div)) +
  geom_violin(size = 0.3, alpha = 0.3, aes(fill = database), colour = NA) + 
  geom_line(aes(group = Sample, alpha=alpha, linewidth = linewidth)) + 
  theme_light() + 
  scale_linewidth_identity() +
  scale_alpha_identity() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values = tool_colours,
                    labels = CCE_names) +
  labs(y = 'Shannon index', x = '', fill = 'Composition estimation')

ggsave('Out/CSHL_poster/Moss_div_MAGs.pdf', bg = 'white', width = 1800, height = 1200, units = 'px', dpi = 180)

#####################
### ALPHA Diversity ##
#######################

dataset <- filter_and_add_samData(
  df = Div_long, 
  ds = 'Moss', 
  Rank = 'Species', 
  ps_list = ps_rare.ls
  ) %>% 
  group_by(database, index) 

Group = 'Compartment'
Databases = c('MPA_db2023','MOTUS','KB20','KB51','SM_genbank-2022.03',"SM_gtdb-rs214-rep")

# Check consistency of testing diversity difference between male/female
test_results <- dataset %>% 
  dplyr::filter(index %in% c('H_1', 'H_2') &
                  database %in% Databases) %>% 
  wilcox_test(as.formula(paste("value ~", Group)), 
              p.adjust.method = NULL) %>% # conservative
  add_significance() %>% 
  dplyr::select(database, index, p.signif)  %>% 
  left_join(dataset, join_by('database', 'index')) %>%
  mutate(p.signif = factor(p.signif, 
                           levels = c('ns','*','**','***', '****')), # reorder factors
         database = factor(database, levels = Databases),
         database = recode(database, !!!CCE_names), #readability
         index = recode(index, !!!Hill_numbers) # ditto
         # value = case_when(index == 'Shannon' ~ log(value),
         #                   index == 'Simpson' ~ 1-1/value,
         #                   TRUE ~ value)
  )

# Plot every database x index combination
test_results %>% # add the p-values to dataset
  ggplot(aes(x = !!sym(Group), y = value, colour = p.signif)) +
  scale_color_discrete() +
  geom_boxplot(linewidth = 0.3, outliers = FALSE) + 
  facet_grid(index~database, scales = 'free_y') +
  #theme_few() +
  theme(
    panel.background = element_rect(fill = "grey95"),  # Sets a lighter grey background
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
    )

ggsave('Out/alpha_P19_Saliva_sex.pdf', bg = 'white',
       width = 3000, height = 1800, units = 'px', dpi = 240)

####################
### BETA Diversity ##
######################
perm <- read_rds('Out/permanova_full.rds')

perm$P19_Gut%>%
  filter(dist == 'Bray-Curtis') %>% 
#  mutate(database = recode(database, !!!CCE_names)) %>% 
  ggplot(aes(x = dist, y = log10(p), colour = database)) +
  geom_beeswarm(aes(size = R2), stroke = 1.5, spacing = 8, shape = 1, 
                method = 'swarm') +  # Use geom_beeswarm to avoid complete overlaps
  facet_grid(. ~ variable, scales = 'free_x') +
  labs(y = expression("log"[10]~"p-value"), 
       colour = 'Tool/database', 
       size = expression("R"^2), 
       shape = 'Distance metric',
       x = '') +
  theme_light() + 
  scale_size_continuous(range = c(0.1,10)) +
  scale_colour_manual(values = tool_colours, labels = CCE_names) +
  theme(axis.text.x = element_blank()) +

  geom_hline(aes(yintercept = log10(0.01), linetype = "p = 0.01"), color = "blue") +
  geom_hline(aes(yintercept = log10(0.05), linetype = "p = 0.05"), color = "red") +
    scale_linetype_manual(name = "", values = c("p = 0.01" = "dashed", "p = 0.05" = "dashed")) +
    guides(
      linetype = guide_legend(order = 1, override.aes = list(color = c("blue", "red"))),
      colour = guide_legend(order = 2),  # Keep the colour legend
      size = guide_legend(order = 4),    # Keep the size legend
      shape = guide_legend(order = 3)    # Keep the shape legend
    )


ggsave('Out/permanova_P19_Gut_species.pdf', bg = 'white', width = 2200, height = 1400, units = 'px', dpi = 180)











