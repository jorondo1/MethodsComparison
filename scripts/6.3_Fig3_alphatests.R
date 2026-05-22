library(pacman)
p_load( magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
        tidyverse, ggh4x, patchwork, ggpubr,
        rstatix)

source("scripts/0_Config.R")
theme_set(theme_light())

#   $$\ $$\         $$$$$$$$\ $$$$$$\  $$$$$$\             $$$$$$\  
#   $$ \$$ \        $$  _____|\_$$  _|$$  __$$\           $$ ___$$\ 
# $$$$$$$$$$\       $$ |        $$ |  $$ /  \__|          \_/   $$ |
# \_$$  $$   |      $$$$$\      $$ |  $$ |$$$$\             $$$$$ / 
# $$$$$$$$$$\       $$  __|     $$ |  $$ |\_$$ |            \___$$\ 
# \_$$  $$  _|      $$ |        $$ |  $$ |  $$ |          $$\   $$ |
#   $$ |$$ |        $$ |      $$$$$$\ \$$$$$$  |$$\       \$$$$$$  |
#   \__|\__|        \__|      \______| \______/ \__|       \______/ 
#                                                                   
##########################
# Alpha diversity tests ###
############################

these_datasets <- c('Moss', 'NAFLD', 'P19_Gut', 'P19_Saliva', 'PD', 'Bee', 'AD_Skin', 'RA_Gut')

these_databases <- c('KB10' ,'KB45', 'KB90', 'SM_RefSeq_20250528', 
                     'KB10_GTDB','KB45_GTDB', 'KB90_GTDB', 
                     'SM_gtdb-rs220-rep', 
                     'MPA_db2022', 'MPA_db2023','MPA_db2025','MOTUS3','MOTUS4')


# Effect size / p-value plot ----------------------------------------------

# Import wilcgox tests, format p
alpha_div_test <- read_rds('Out/_Rdata/alpha_div.RDS')[['wilcox_tests']] %>% # conservative
  mutate(p.signif = case_when(
    p < 0.001 ~ 'p < 0.001',
    p < 0.01 ~ 'p < 0.01',
    p < 0.05 ~ 'p < 0.05',
    TRUE ~ 'p ≥ 0.05'),
    n = n1+n2) %>% 
  dplyr::select(Dataset, Database, Index,effsize, magnitude, p, p.signif, n) %>% 
  mutate(
    p.signif = factor(
      p.signif, 
      levels = c('p ≥ 0.05', 'p < 0.05', 'p < 0.01', 'p < 0.001' )))

# dataset factor levels with n inside
dataset_names_n_tibble <- alpha_div_test %>% 
  distinct(Dataset, n) %>% 
  group_by(Dataset) %>% 
  filter(n==max(n)) %>% 
  left_join(tibble(
    Dataset = names(dataset_names),
    Dataset_name = unname(dataset_names)
  ), by = "Dataset") %>% 
  mutate(Dataset_name_n = paste0(Dataset_name, "\n(n=",n,")" ))

# values:
dataset_names_n_lookup <- dataset_names_n_tibble$Dataset_name_n
# add standard names:
names(dataset_names_n_lookup) <- dataset_names_n_tibble$Dataset 
# reorder and fluhs NAs
dataset_names_n_lookup %<>% .[names(dataset_names)] %>% .[!is.na(.)]

which_indices <- c('Richness', 'Shannon', 'Hill_2', 'Tail')

# plot data
wilcox.pdat <- CCE_metadata %>% 
  dplyr::select(Database, CCE_approach, short_alpha_2, short_alpha_3) %>% 
  right_join(alpha_div_test, by = 'Database') %>% 
  filter(Dataset %in% these_datasets, 
         Database %in% these_databases,
         Index %in% which_indices) %>% 
  mutate(Dataset = factor(recode(Dataset, !!!dataset_names_n_lookup), levels = dataset_names_n_lookup),
         CCE_approach = factor(CCE_approach, levels = c("D2M", "D2D"))) %>% 
  split(., .['Index']) # split by group 

wilcox.plots <- imap(wilcox.pdat, function(wilcox.df, Index_name) {
  
  plot <- wilcox.df %>% 
    ggplot(aes(x = Database, y = Dataset, size = effsize, colour = p)) +
    geom_point() +
    facet_nested(Dataset ~ CCE_approach + short_alpha_2 + short_alpha_3, scales = 'free') +
    scale_colour_steps2(
      low = "#3d0173", mid = "#f0f0f0", high = "#cccccc",
      midpoint = log10(0.05), 
      trans = "log10",
      breaks = c(0.001, 0.01, 0.05, 0.1, 0.5, 1),
      limits = c(0.001, 1), 
      oob = scales::squish
    ) +
    scale_size_area(max_size = 12) +
    theme_light()+
    labs(size = "Rank-biserial correlation", colour = 'p-value') +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = 'bottom',
      legend.title.position = 'top',
      legend.key.width = unit(1, "cm"),
      !!!strip_theme) 

})

imap(wilcox.plots, function(plot, idx){
  ggsave(paste0('Out/Manuscript/3.1.alpha_wilcox_',idx, '.pdf'),
         plot = plot,
         bg = 'white', width = 2500, height = 2000,
         units = 'px', dpi = 280)
})

# Tables ------------------------------------------------------------------

imap(wilcox.pdat, function(wilcox.df, Index_name) {
  
  wilcox.df %>% 
    dplyr::select(Database, Dataset, effsize, p)%>% 
    mutate(value = paste0(format(effsize, digits=1), 
                          " [p=", round(p, digits = 4), "]"),
           .keep = 'unused') %>% 
    pivot_wider(names_from = 'Database', values_from = value) %>% 
    kable() %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>% 
    save_kable(paste0('Out/Manuscript/wilcox_',Index_name, '.html'))
})

# Violin plots  -----------------------------------------------------------

div.data <- read_rds('Out/_Rdata/alpha_div.RDS')[['plot_data']] %>% 
  filter(Dataset == "P19_Gut", 
         Database %in% c("MPA_db2025", "KB45"),
         Index == 'Shannon') %>% 
  mutate(Database = recode(Database, !!!CCE_names)) %>% 
  mutate(Database = str_replace(Database, ' \\(', '\n('))

stat.test <- div.data %>%
  group_by(Database) %>% 
  wilcox_test(Index_value~ Grouping_var) %>%
  add_xy_position(x = "Grouping_var") %>% 
  mutate(p_equals = paste0("p=", round(p,5)),
         Index = "Shannon index")

Shannon_violin <- div.data %>% 
  mutate(Diarrhea = ifelse(Grouping_var=="A", 'Yes', 'No'),
         Index = "Shannon index") %>% 
  # plot:
  ggplot(aes(x = Diarrhea, y = Index_value, fill = Diarrhea)) + 
  geom_violin(draw_quantiles = 0.5)+
  geom_jitter(width = 0.2) +
  stat_pvalue_manual(stat.test, label = "p_equals") +
  facet_nested(~Index+Database) +
  labs(y = 'Shannon index') +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        legend.title.position = 'top',
        !!!strip_theme
  )

# --COWPLOT:  Label each plot
p_A <- ggdraw() +
  draw_plot(wilcox.plots$Shannon + 
              theme(legend.position = 'bottom',
                    legend.title.position = 'top',
                    legend.key.width = unit(1, "cm"))) +
  draw_label("A.", x = 0.02, y = 0.985, hjust = 0, vjust = 1,
             fontface = "bold", size = 14)

p_B <- ggdraw() +
  draw_plot(Shannon_violin) +
  draw_label("B.", x = 0.0, y = 0.985, hjust = 0, vjust = 1,
             fontface = "bold", size = 14)

# Assemble — AAAAABB = 5 units A, 2 units B
plot_grid(
  p_A, p_B,
  ncol = 2,
  rel_widths = c(7, 2)
  )

ggsave(paste0('Out/Manuscript/3.2.alpha_wilcox_example.pdf'),
       bg = 'white', width = 2500, height = 2000,
       units = 'px', dpi = 260)
# 
# # Diversity indices
# alpha_div <- read_rds('Out/_Rdata/alpha_div.RDS')[['plot_data']] %>% 
#   filter(Dataset %in% these_datasets, Database %in% these_databases) %>% 
#   left_join(CCE_metadata %>% select(Database, CCE_approach, short_alpha_2, short_alpha_3), by = 'Database') %>% 
#   # add p values
#   left_join(alpha_div_test, by = c('Dataset', 'Database', 'Index')) %>% 
#   mutate(
#     Database = factor(Database, levels = these_databases),
#     label = paste0(
#       'p=',
#       ifelse(
#         p < 0.01, 
#         format(p, scientific = TRUE, digits = 2), 
#         round(p, 3))
#       )
#     )
# 
# 
# # Plot and save :
# make_plots <- function(desc, idx) {
#   
#   # Skip KB45 for Richness
#   database_subset <- if(idx %in% c('Richness', 'Tail') ){
#     these_databases[which(!these_databases %in% c("KB45", "KB10", 'KB10_GTDB', 'KB45_GTDB'))]
#   } else {these_databases}
#   
#   # Create one-row-per-facet label data
#   alpha_div_labels <- alpha_div %>%
#     filter(Index == idx, Database %in% database_subset) %>%
#     distinct(Dataset, Database, CCE_approach, short_alpha_2, short_alpha_3, label, p.signif)
#   
#   alpha_div %>% 
#     filter(Index == idx, Database %in% database_subset) %>%
#     ggplot(aes(x = Grouping_var, y = Index_value, fill = p.signif)) +
#     geom_violin(linewidth =0.2, draw_quantiles = c(0.50)) +
#     facet_nested(Dataset ~ CCE_approach + short_alpha_2 + short_alpha_3, scales = 'free') +
#     geom_text(
#       data = alpha_div_labels,
#       aes(x = 0.5, y = Inf, label = label),
#       hjust = 0, vjust = 1.5,
#       size = 3, color = "black",
#       inherit.aes = FALSE
#     )+
#     labs(fill = 'p-value', x = 'Group', y = desc) +
#     scale_fill_manual(values = c(
#       "#b7cdd0",
#       "#3ba7e5",
#       "#8d74ce",
#       "#83ab56")) +
#     theme(
#       legend.position = c(0.9,0.75),
#       legend.title = element_blank(),
#       strip.background = element_rect(fill = 'grey50'),
#       axis.text = element_text(size = 8),
#       legend.background = element_rect(
#         fill = "white",        
#         color = "black",       # border
#         linewidth = 0.2        # Border
#       )) 
# }
# 
# axis_desc_tests <- c(
#   Richness = 'Species richness',
#   Shannon = "Shannon's diversity index",
#   Simpson = "Simpson's dominance index",
#   Tail = "Tau statistic"
# )
# 
# plots <- imap(axis_desc_tests, make_plots)
# 
# imap(plots, function(p, idx) {
#   ggsave(plot = p, paste0('Out/Manuscrit/alpha_',idx, '_tests.pdf'),
#          bg = 'white', width = 2000, height = 2400,
#          units = 'px', dpi = 220)
# })
# 
# # Power, p-values and n for each method
# 
# alpha_div
