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

dataset_names_n_lookup <- read_rds('Out/_Rdata/dataset_names_n_lookup_alphadiv.RDS')

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
      !!!strip_theme)  +
    guides(
      size = guide_legend(order = 1)  # appears second = rightmost
    )
  

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
        axis.ticks.x = element_blank(),
        legend.position = 'bottom',
        legend.text.position = 'bottom',
        
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

# Assemble 
plot_grid(
  p_A, p_B,
  ncol = 2,
  rel_widths = c(7, 2)
  )

ggsave(paste0('Out/Manuscript/3.2.alpha_wilcox_example.pdf'),
       bg = 'white', width = 2500, height = 2000,
       units = 'px', dpi = 260)

ggsave('Out/Manuscript/PLOS/Fig4.eps',
       bg = 'white', device = 'eps',
       width = 2900, height = 2320, units = 'px')

