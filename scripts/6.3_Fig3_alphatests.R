library(pacman)
p_load( magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
        tidyverse, ggh4x,
        rstatix)

ps_rare.ls <- read_rds('Out/_Rdata/ps_rare.ls.rds')
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

# Test between groups
# Simple violins with pvalues
these_datasets <- c('Moss', 'NAFLD', 'P19_Gut', 'P19_Saliva', 'PD', 'Bee', 'AD_Skin')

these_databases <- c('KB45', 'KB90', 'SM_RefSeq_20250528', 
                     'KB45_GTDB', 'KB90_GTDB', 'SM_gtdb-rs220-rep',
                     'MPA_db2023','MOTUS')

alpha_div <- read_rds('Out/_Rdata/alpha_div.RDS')[['plot_data']] %>% 
  filter(Dataset %in% these_datasets) %>% 
  mutate(Database = factor(Database, levels = names(tooldb_colours))) %>% 
  left_join(CCE_metadata, by = 'Database')

axis_desc_tests <- c(
  Richness = 'Species richness',
  Shannon = "Shannon's diversity index",
  Simpson = "Simpson's dominance index"
)

# Prepare data
alpha_div_test <- read_rds('Out/_Rdata/alpha_div.RDS')[['wilcox_tests']] %>% # conservative
  mutate(p.signif = case_when(
    p < 0.001 ~ 'p < 0.001',
    p < 0.01 ~ 'p < 0.01',
    p < 0.05 ~ 'p < 0.05',
    TRUE ~ 'p ≥ 0.05')) %>% 
  select(Dataset, Database, Index, p, p.signif) %>% 
  mutate(
    p.signif = factor(
      p.signif, 
      levels = c('p ≥ 0.05', 'p < 0.05', 'p < 0.01', 'p < 0.001' )),
    Methode = case_when(
      str_detect(Database, 'KB') | str_detect(Database, 'SM') ~ 'DNA-to-DNA',
      str_detect(Database, 'MPA') | str_detect(Database, 'MOTUS') ~ 'DNA-to-marker'),
    RefDB = case_when(
      str_detect(Database, regex("GTDB", ignore_case = TRUE)) ~'GTDB',
      Methode == 'DNA-to-marker'  ~ ' ',
      TRUE ~ "RefSeq"))

alpha_labels <- function(test_df, IDX) {
  filter(test_df,
         Index == IDX 
         & Database %in% these_databases
         & Dataset %in% these_datasets) %>% 
    mutate(
      Database = factor(Database, levels = these_databases),
      label = paste0('p=',
                     ifelse(p < 0.01, 
                            format(p, scientific = TRUE, digits = 2), 
                            round(p, 3)))
    )
}

# Plot and save :
plots <- imap(axis_desc_tests, function(desc, idx) {
  
  # Skip KB45 for Richness
  database_subset <- if(idx == 'Richness') {
    these_databases[!str_detect(these_databases, 'KB45')]
  } else {these_databases}
  
  # Set p value labels
  alpha_labs <- alpha_div_test %>% 
    filter(Database %in% database_subset) %>% 
    alpha_labels(., idx)
  
  # Plot
  alpha_div %>% 
    filter(Index == idx 
           & Database %in% database_subset
    ) %>% 
    left_join(alpha_labs, 
              by = c('Dataset', 'Database', 'Index')) %>% 
    mutate(Database = factor(Database, levels = database_subset)) %>% 
    ggplot(aes(x = Grouping_var, y = Index_value, fill = p.signif)) +
    geom_violin(linewidth =0.2, draw_quantiles = c(0.50)) +
    geom_text(
      data = alpha_labs,
      aes(x = 0.5, y = Inf, label = label, 
          hjust = 0,
          vjust = 14), # controls (weirdly) the y position of text label
      size = 3, color = "black"
    ) + 
    facet_nested(
      Dataset ~ Methode + RefDB + Database, scales = 'free',
      # Relabel facets :
      labeller = labeller(Database = as_labeller(
        setNames(CCE_metadata$MethodNameParam, CCE_metadata$Database)
      ))) +
    labs(fill = 'p-value', x = 'Group', y = desc) +
    scale_fill_manual(values = c(
      "#b7cdd0",
      "#3ba7e5",
      "#8d74ce",
      "#83ab56")) +
    theme(
      legend.position = c(0.9,0.75),
      legend.title = element_blank(),
      strip.background = element_rect(fill = 'grey50'),
      axis.text = element_text(size = 8),
      legend.background = element_rect(
        fill = "white",        
        color = "black",       # border
        linewidth = 0.2        # Border
      )) 
})

imap(plots, function(p, idx) {
  ggsave(plot = p, paste0('Out/memoire/alpha_',idx, '_tests.pdf'),
         bg = 'white', width = 2000, height = 2400,
         units = 'px', dpi = 220)
})

# Power, p-values and n for each method
