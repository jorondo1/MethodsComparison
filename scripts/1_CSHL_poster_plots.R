library(pacman)
p_load(
  # syntax
  magrittr, tidyverse, purrr, furrr, foreach, doParallel, parallel,
  # metrics and stats
  phyloseq, DESeq2, vegan, rstatix,
  # plotting :
  ggbeeswarm2, patchwork, grid, ggh4x, ggthemes, ComplexHeatmap, UpSetR)

ps.ls <- read_rds('Out/ps_full.ls.rds')
ps_rare.ls <- read_rds('Out/ps_rare.ls.rds')

#####################
### ALPHA Diversity ##
#######################

dataset <- filter_and_add_samData(
  df = Div_long, 
  ds = 'P19_Saliva', 
  Rank = 'Species', 
  ps_list = ps_rare.ls
  ) %>% 
  group_by(database, index) 

Group = 'sex'
Databases = c('MPA_db2022','MPA_db2023','MOTUS','KB20','KB51',"SM_gtdb-rs214-rep")

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
         database = recode(database, !!!CCE_names), #readability
         index = recode(index, !!!Hill_numbers), # ditto
         # value = case_when(index == 'Shannon' ~ log(value),
         #                   index == 'Simpson' ~ 1-1/value,
         #                   TRUE ~ value)
  ) %>% mutate(database = factor(database, levels = CCE_names))

# Plot every database x index combination
test_results %>% # add the p-values to dataset
  ggplot(aes(x = !!sym(Group), y = value, colour = p.signif)) +
  scale_color_discrete() +
  geom_boxplot(linewidth = 0.3, outliers = FALSE) + 
  facet_grid(index~database, scales = 'free_y') +
  theme_light() +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = 'bottom',
    legend.text.position = 'bottom',
    legend.title.position = 'top',
    legend.title = element_text(hjust = 0.5),
    strip.text = element_text(color = "black")
    ) +
  scale_colour_manual(values = c('red3', 'blue3', 'green3')) +
  labs(y = 'Effective species (Hill numbers)', x = '', colour = 'Wilcoxon')

ggsave('Out/CSHL_poster/alpha_P19_Saliva_sex.pdf', bg = 'white',
       width = 2500, height = 1800, units = 'px', dpi = 260)

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



################
### Mosses ######
##################

Div_long <- read_rds('Out/Diversity_long_full.rds') %>% 
  mutate(index_div = case_when(index == 'H_1' ~ log(value),
                               index == 'H_2' ~ 1-1/value, 
                               TRUE ~ value))

## Maybe one comparative diversity plot using lines only?
# Find most variable lines per group
Div_long_filt <- Div_long %>% 
  filter(dataset %in% c('Moss') & 
           index == 'H_1' &
           Rank == 'Species' &
           !database %in% c('KB51')) 

# set.seed(38); most_variable <- Div_long_filt %>% 
#   dplyr::select(index, Rank, dataset, Sample) %>%
#   unique %>% 
#   mutate(random = runif(n(), min = 0, max = 1)) %>% 
#   filter(random > 0.94) %>% 
#   mutate(alpha = 1, linewidth = 0.3) 

most_variable <- Div_long_filt %>% 
  dplyr::filter(!database %in% c('SM_gtdb-rs214-full', 'KB20')) %>% 
  group_by(Sample) %>% 
  dplyr::summarise(var_index = var(value)) %>% 
  arrange(var_index) %>% tail(n=12) %>% 
  mutate(alpha = 1, linewidth = 0.3) 
  
# The most variable lines could be thicker...?
Div_long_filt %>% 
  left_join(most_variable, by='Sample') %>% 
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


###########
#### DAA ###
#############
DAA_files <- 'Out/DAA*/*.tsv'
DAA <- Sys.glob(DAA_files) %>% map(read_tsv) %>% list_rbind

# NAFLD subset of significant taxa

upset_by_db <- function(ds, db) {
  wide_DAA <- ds %>% 
    filter(database == db) %>% 
  ### Wide matrix for upset plot 
    dplyr::select(Taxon, DAA_tool) %>% 
    distinct() %>% 
    mutate(present = 1) %>% 
    pivot_wider(names_from = DAA_tool, values_from = present, values_fill = 0) 
  
  # Simple upset plot 
  upset(wide_DAA %>% as.data.frame, 
        sets = names(wide_DAA)[-1],
        order.by = "freq",
        mainbar.y.max = 65,
        nintersects = 14,
        set_size.scale_max = 160,
        mb.ratio = c(0.5,0.5))
  
}
Moss_Species_sig <- DAA %>% 
  filter(taxRank == 'Species' &
           dataset == 'Moss' &
           DAA_tool %in% c('ANCOMBC2', 'Aldex2', 'ZicoSeq', 'corncob', 'radEmu') &
           adj.p <0.05) 

upset_by_db(Moss_Species_sig, 'SM_gtdb-rs214-rep')
upset_by_db(Moss_Species_sig, 'SM_gtdb-rs214-rep_MAGs')

Moss_Species_sig %>% 
  group_by(Taxon, database) %>% 
  dplyr::summarise(n = n(), .groups='drop') %>% 
  filter(n>3) %>% 
  group_by(database) %>% 
  dplyr::summarise(n = n())
  
