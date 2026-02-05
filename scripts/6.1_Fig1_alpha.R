library(pacman)
p_load( magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
        tidyverse, kableExtra, gghalves,
        rstatix)

ps_rare.ls <- read_rds('Out/_Rdata/ps_rare.ls.rds')
source("scripts/0_Config.R")
theme_set(theme_light())

#   $$\ $$\         $$$$$$$$\ $$$$$$\  $$$$$$\              $$\   
#   $$ \$$ \        $$  _____|\_$$  _|$$  __$$\           $$$$ |  
# $$$$$$$$$$\       $$ |        $$ |  $$ /  \__|          \_$$ |  
# \_$$  $$   |      $$$$$\      $$ |  $$ |$$$$\             $$ |  
# $$$$$$$$$$\       $$  __|     $$ |  $$ |\_$$ |            $$ |  
# \_$$  $$  _|      $$ |        $$ |  $$ |  $$ |            $$ |  
#   $$ |$$ |        $$ |      $$$$$$\ \$$$$$$  |$$\       $$$$$$\ 
#   \__|\__|        \__|      \______| \______/ \__|      \______|

####################
# Alpha diversity ###
######################

# Reorder factors, filter for datasets
these_datasets <- c('Moss', 'NAFLD', 'P19_Gut', 'P19_Saliva', 'PD', 'Bee', 'AD_Skin')

alpha_div <- read_rds('Out/_Rdata/alpha_div.RDS')[['plot_data']] %>% 
  filter(Dataset %in% these_datasets) %>% 
  mutate(Database = factor(Database, levels = names(tooldb_colours))) %>% 
  left_join(CCE_metadata, by = 'Database')

# PLOT hill_1 variation for methods most equivalent in terms of number of species
these_databases <- c('MPA_db2023', 'MOTUS', 
                     'KB45', 'KB45_GTDB', 
                     'SM_gtdb-rs220-rep', 'SM_RefSeq_20250528')

axis_desc <- c(
  Richness = 'Hill number of order 0 (Number of species)',
  Hill_1 = 'Hill number of order 1 (effective number of equally abundant species)',
  Hill_2 = 'Hill number of order 2 (effective number of dominant species)'
)

div_comparison.pdat <- alpha_div %>% 
  filter(Database %in% these_databases,
         Dataset %in% these_datasets) %>% 
  mutate(Database = factor(Database, levels = these_databases))

# Keep only paired samples for each method comparison
keep_paired_samples <- function(df, idx) {
  require(magrittr)
  
  df %<>% filter(Index == idx)
  
  # Keep only samples who make it across both methods
  sample_db_sets <- df %>% 
    group_by(Taxonomy, Sample) %>% 
    summarise(n = n(), .groups = 'drop') %>% 
    filter(n == 2) %>% 
    select(-n)
  
  sample_subset <- left_join(
    x = sample_db_sets,
    y = df,
    by = c('Taxonomy', 'Sample')
  )
  
  # Count samples
  counts <- sample_subset %>% 
    group_by(Taxonomy) %>% 
    summarise(n = n()/2) %>% 
    deframe() # make it a named vector
  
  # Add sample count by facet
  sample_subset %>% 
    mutate(
      Facet = case_when(
        Taxonomy == 'Tool-specific' ~ paste0('A. DNA-to-marker methods (n = ', counts['Tool-specific'], ')'),
        Taxonomy == 'GTDB' ~ paste0('B. GTDB 220 (n = ', counts['GTDB'], ')'),
        Taxonomy == 'NCBI' ~ paste0('C. RefSeq 2024-12-28 (n = ', counts['NCBI'], ')'))
    )
}


# PANEL 1 : distribution of index values ----------------------------------

alpha_plots <- imap(axis_desc, function(desc, idx) {
  
  dat <- div_comparison.pdat %>% 
    keep_paired_samples(idx = idx) %>% 
    filter(Index == idx)
  
  dat %>% 
    ggplot(aes(x = Database, y = Index_value)) +
    geom_half_violin(
      data = . %>% filter(str_detect(Database, 'KB') | str_detect(Database, 'MPA')),
      aes(fill = Tool), 
      linewidth = 0.3,
      side = 'l',
      draw_quantiles = 0.5) + 
    geom_half_violin(
      data = . %>% filter(str_detect(Database, 'SM') | str_detect(Database, 'MOTUS')),
      aes(fill = Tool), 
      linewidth = 0.3,
      side = 'r',
      draw_quantiles = 0.5) + 
    geom_line(aes(group = Sample), alpha = 0.5, linewidth = 0.08) +
    facet_wrap(~Facet, scales = 'free') +
    scale_fill_manual(values = tool_colours, breaks = c('MetaPhlAn4', 'mOTUs3', 'Kraken2+Bracken', 'Sourmash gather')) +
    theme_light() +
    theme(
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = 'bottom',
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.box.spacing = unit(-0.5, "lines"),
      strip.background = element_rect(fill = 'grey50'),
      strip.text.x.top = element_text(
        angle = 0, hjust = 0, size = 12),
      legend.background = element_rect(
        fill = "white",        # White background
        color = "black",       # Black border
        linewidth = 0.3        # Border thickness
      )
    ) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(fill = 'Tool', x = '', y = desc)
  
})


# SAVE PLOTS
imap(alpha_plots, function(plot, idx) {
  
  ggsave(plot = plot, 
         paste0('Out/memoire/alpha_',idx, '_comparison.pdf'),
         bg = 'white', width = 2200, height = 1400,
         units = 'px', dpi = 220)
  
  ggsave(plot = plot, 
         paste0('Out/ISMB2025/alpha_',idx, '_comparison.pdf'),
         bg = 'white', width = 2300, height = 1200,
         units = 'px', dpi = 230)
  
})


# PANEL 2 : distribution of differences -----------------------------------

centered_differences <- div_comparison.pdat %>% 
  keep_paired_samples(idx = "Hill_1") %>% 
  filter(Index == "Hill_1") %>% 
  #  group_by(Facet) %>% 
  #  mutate(Centered_value = Index_value - median(Index_value)) %>% 
  group_by(Taxonomy, Dataset, Facet, Sample) %>%
  #  summarise(differences = last(Centered_value)-first(Centered_value),
  summarise(differences = (last(Index_value)-first(Index_value))/last(Index_value),
            .groups = 'drop') %>%
  group_by(Facet) %>% 
  mutate(centered_diffs = differences - median(differences)
  ) %>% 
  ungroup()

centered_differences %>% 
  #  filter(centered_diffs<2) %>% 
  
  ggplot(aes(x = centered_diffs, fill = Facet)) +
  geom_density(alpha = 0.4) +
  labs(x = "Median-centered changes in Hill diversity (order 1)") +
  theme(
    legend.position = 'bottom'
  )

# Variances
centered_differences %>% 
  group_by(Facet) %>% 
  summarise(var_cdiff = var(centered_diffs)) 

centered_differences %>% 
  group_by(Facet) %>% 
  mutate(squared_dev = (differences - median(differences))^2) %>% 
  ungroup() %>% 
  arrange(Sample) %>% 
  wilcox_test(squared_dev~Facet)


# Signed rank test on square deviations from the median
# suggested by https://claude.ai/chat/55cf7920-1402-4e40-85ba-074270266e55
centered_differences %>% 
  # only keep samples evaluated by all 3 pairs
  group_by(Sample) %>% 
  #filter(n()==3) %>%
  filter(Taxonomy != 'Tool-specific') %>% 
  filter(n()==2) %>% 
  group_by(Facet) %>% 
  mutate(squared_dev = (differences - median(differences))^2) %>% 
  ungroup() %>% 
  arrange(Sample) %>% 
  wilcox_test(squared_dev~Facet, paired = TRUE)



# Sample diversity variations tables --------------------------------------

quantify_div_variation <- function(df, ds1, ds2, idx) {
  dir_change <- df %>% 
    select(Sample, Database, Index_value, Dataset, Index) %>% 
    filter(Database %in% c(ds1,ds2)
           & Index == idx) %>% 
    pivot_wider(names_from = Database,
                values_from = Index_value) %>% 
    mutate(change = .[[ds1]] - .[[ds2]], # Difference between 1st and 2nd tool
           higher_with = case_when(
             sign(change) == 1 ~ ds1,
             sign(change) == -1 ~ ds2,
             TRUE ~ 'None')
    )%>%  # Dynamic column name in mutate()
    filter(!is.na(change))
  
  # Summarise both by dataset and overall
  summary.ls <- list(
    dir_change %>% group_by(Dataset, higher_with),
    dir_change %>% group_by(higher_with)
  ) %>% 
    lapply(function(df) {
      summarise(
        df,
        median_change = abs(median(change)),
        mad_c = mad(change),
        min_c = min(abs(change)),
        max_c = max(abs(change)),
        rcv_c = mad_c/median_change,
        n = n(), .groups = 'drop')
    }
    ) 
  summary.ls[[2]] %<>% mutate(Dataset = 'All datasets')
  
  list_rbind(summary.ls) %>% 
    # Compute percentage samples
    group_by(Dataset) %>% 
    mutate(prop = 100*n / sum(n),
           across(where(is.numeric), ~ round(.x, 2))) %>% 
    kable(align = "l", 
          caption = paste0(" ", idx, ": ",ds1, ' values minus ', ds2, ' values')) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) 
}

for (idx in c('Richness', 'Hill_1', 'Hill_2')) {
  quantify_div_variation(
    alpha_div, idx = idx,
    'KB10_GTDB', 'SM_gtdb-rs220-rep') %>% 
    save_kable(paste0('Out/memoire/tables/alpha_', idx,'_GTDB10.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'KB45_GTDB', 'SM_gtdb-rs220-rep') %>% 
    save_kable(paste0('Out/memoire/tables/alpha_', idx,'_GTDB45.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'KB90_GTDB', 'SM_gtdb-rs220-rep') %>% 
    save_kable(paste0('Out/memoire/tables/alpha_', idx,'_GTDB90.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'KB10', 'SM_RefSeq_20250528') %>% 
    save_kable(paste0('Out/memoire/tables/alpha_', idx,'_RefSeq10.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'KB90', 'SM_RefSeq_20250528') %>% 
    save_kable(paste0('Out/memoire/tables/alpha_', idx,'_RefSeq90.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'KB45', 'SM_RefSeq_20250528') %>% 
    save_kable(paste0('Out/memoire/tables/alpha_', idx,'_RefSeq45.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'MPA_db2023', 'MOTUS') %>% 
    save_kable(paste0('Out/memoire/tables/alpha_', idx,'_Markers.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'SM_gtdb-rs220-rep', 'SM_RefSeq_20250528') %>% 
    save_kable(paste0('Out/memoire/tables/alpha_', idx,'_SM.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'KB45_GTDB', 'KB45') %>% 
    save_kable(paste0('Out/memoire/tables/alpha_', idx,'_KB.html'))
}

# Mean Dataset alphadiv range across methods
alpha_div %>% group_by(Database, Dataset, Taxonomy) %>% 
  filter(Index %in% c('Hill_1')
         & Database %in% these_databases) %>%
  # Mean tool div by dataset
  summarise(mean = mean(Index_value), .groups = 'drop') %>% 
  # Mean 
  group_by(Dataset, Taxonomy) %>% 
  summarise(min = min(mean), max = max(mean), .groups = 'drop') %>% 
  group_by(Taxonomy) %>% 
  mutate(fold_increase = max / min) %>% 
  summarise(mean_fold = mean(fold_increase),
            sd_fold = sd(fold_increase))
