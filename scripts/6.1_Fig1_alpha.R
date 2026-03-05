library(pacman)
p_load( magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
        tidyverse, kableExtra, gghalves,
        rstatix)

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

# TEST_STARTPOINT for dropout tests

####################
# Alpha diversity ###
######################

# Reorder factors, filter for datasets
these_datasets <- c('Moss', 'NAFLD', 'P19_Gut', 'P19_Saliva', 'PD', 'Bee', 'AD_Skin', 'RA_Gut')

alpha_div <- read_rds('Out/_Rdata/alpha_div.RDS')[['plot_data']] %>% 
  filter(Dataset %in% these_datasets) %>% 
  mutate(Database = factor(Database, levels = names(tooldb_colours))) %>% 
  left_join(CCE_metadata, by = 'Database')

# PLOT hill_1 variation for methods most equivalent in terms of number of species
these_databases <- c('MPA_db2025', 'MOTUS4', 
                     'KB45', 'KB45_GTDB', 
                     'SM_gtdb-rs220-rep', 'SM_RefSeq_20250528')

axis_desc <- c(
  Richness = 'Hill number of order 0 (Number of species)',
  Shannon = 'Shannon entropy',
  Hill_1 = 'Hill number of order 1 (effective number of equally abundant species)',
  Simpson = 'Simpson index',
  Hill_2 = 'Hill number of order 2 (effective number of dominant species)',
  Tail = 'Tail statistic (rank-based diversity)'
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
  
  sample_subset %>% 
    mutate(
      paired_N = counts[Taxonomy]
    )
  
}

alpha_paired <- map(names(axis_desc), function(idx) {
  div_comparison.pdat %>% 
    keep_paired_samples(idx = idx)
}) %>% setNames(names(axis_desc))


# PANEL 1 : distribution of index values ----------------------------------

alpha_distr.plots <- imap(alpha_paired, function(dat, idx) {
  
  dat %>% 
    mutate(
      Facet = case_when(
        Taxonomy == 'Tool-specific' ~ paste0('A. DNA-to-marker methods (n = ', paired_N, ')'),
        Taxonomy == 'GTDB' ~ paste0('B. GTDB 220 (n = ', paired_N, ')'),
        Taxonomy == 'NCBI' ~ paste0('C. RefSeq 2024-12-28 (n = ', paired_N, ')'))
    ) %>% 
    # Plot :
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
    scale_fill_manual(values = tool_colours, breaks = c('MetaPhlAn4', 'mOTUs4', 'Kraken2+Bracken', 'Sourmash gather')) +
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
        fill = "white",    
        color = "black",   
        linewidth = 0.3    
      )
    ) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(fill = 'Tool', x = '', y = axis_desc[idx])
  
})


# SAVE PLOTS
imap(alpha_distr.plots, function(plot, idx) {
  
  ggsave(plot = plot, 
         paste0('Out/Manuscript/alpha_',idx, '_comparison.pdf'),
         bg = 'white', width = 2200, height = 1400,
         units = 'px', dpi = 210)
  # 
  # ggsave(plot = plot, 
  #        paste0('Out/ISMB2025/alpha_',idx, '_comparison.pdf'),
  #        bg = 'white', width = 2300, height = 1200,
  #        units = 'px', dpi = 230)
  
})

# PANEL 2 : distribution of differences -----------------------------------
# Intéressant : encore plus variable/prononcé avec Hill_2

alpha_diff.pdat <- map(alpha_paired, function(dat) {
  
  dat %>% 
    # Text for legend titles presented differently:
    mutate(
      Facet = case_when(
        Taxonomy == 'Tool-specific' ~ paste0('MetaPhlAn4 (2025) vs. mOTUs4 (n = ', paired_N, ')'),
        Taxonomy == 'GTDB' ~ paste0('Kraken vs. Sourmash (GTDB 220; n = ', paired_N, ')'),
        Taxonomy == 'NCBI' ~ paste0('Kraken vs. Sourmash (RefSeq 2024-12-28; n = ', paired_N, ')'))
    ) %>% 
    select(Taxonomy, Dataset, Database, Facet, Sample, Index_value) %>% 
    
    # Scale indices (median/mad) within dataset/database combo:
    group_by(Dataset, Database) %>% 
    mutate(scaled_index = (Index_value-median(Index_value))/mad(Index_value)) %>% 
    ungroup() %>% 
    
    # Pairwise differences in scaled diversity :
    group_by(Taxonomy, Dataset, Facet, Sample) %>%
    
    # each group will have 2 by design (see these_databases vector)
    # so we can subtract the first one from the last one within group
    # if they are ordered :
    arrange(Taxonomy, Dataset, Facet, Sample) %>% 
    summarise(differences = first(scaled_index)-last(scaled_index),
              .groups = 'drop') 
}) 

distr_alpha_diffs.plots <- imap(alpha_diff.pdat, function(plot.dat, plot_name){
  
  plot.dat %>% 
    # Format facet labels with sample counts 
    #  filter(centered_diffs<2) %>% 
    ggplot(aes(x = differences, fill = Facet)) +
    geom_density(alpha = 0.4, linewidth = 0.2) +
    theme(
      legend.position = 'bottom'
    ) +
    labs(subtitle = paste0(LETTERS[which(names(alpha_diff.pdat) == plot_name)], ": ", plot_name),
         fill = "Tool pair comparison") +
    theme(plot.caption = element_text(hjust = 0),
          axis.title = element_blank()) 
  
})

# Combine in 2 columns
wrap_plots(distr_alpha_diffs.plots, ncol = 2) +
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom',
        legend.title = element_blank()) &
  plot_annotation(title = "Density of changes in scaled diversity indices")

ggsave(paste0('Out/Manuscript/alpha_diff_distribution.pdf'),
       bg = 'white', width = 2500, height = 1400,
       units = 'px', dpi = 260)
# 
# # Variances
# centered_differences %>% 
#   group_by(Facet) %>% 
#   summarise(var_cdiff = var(differences)) 
# 
# # Do variances differ significantly?
# centered_differences %>% 
#   group_by(Facet) %>% 
#   mutate(squared_dev = (differences - median(differences))^2) %>% 
#   ungroup() %>% 
#   arrange(Sample) %>% 
#   wilcox_test(squared_dev~Facet)

# Spearman correlation of ranks -------------------------------------------

spearman.raw <- imap(alpha_paired, function(dat, idx) {
  
  dat %>% 
    select(Taxonomy, Database, Sample, Index_value) %>% 
    group_split(Taxonomy) %>%
    map_dfr(function(x) {
      x %>%
        pivot_wider(names_from = Database, values_from = Index_value) %>%
        select(-Sample, -Taxonomy) %>% 
        cor_test(vars = everything(), method = "spearman")  %>%
        transmute(
          Taxonomy = unique(x$Taxonomy),
          Index = idx,
          cor = cor,
          p = p)
    })
}) %>% list_rbind() 

spearman.table <- spearman.raw %>% 
  # hill numbers are monotonic transformations, hence redundant
  filter(!Index %in% c("Hill_1", "Hill_2")) %>% 
  mutate(
    Approach = case_when(
      Taxonomy == 'Tool-specific' ~ paste0('MetaPhlAn4 (2025) vs. mOTUs4'),
      Taxonomy == 'GTDB' ~ paste0('Kraken vs. Sourmash (GTDB 220)'),
      Taxonomy == 'NCBI' ~ paste0('Kraken vs. Sourmash (RefSeq 2024-12-28)')),
    .keep = 'unused', .after = Index
  ) %>% 
  pivot_wider(names_from = Index, values_from = cor, id_cols = Approach) %>% 
  kable() #%>% 
  #kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

#save_kable(spearman.table, 'Out/Manuscript/alpha_corr.html')

# Drop-out test PD --------------------------------------------------------

# >>>> TEST_BREAKPOINT # used to subset the script (stops here)
# save current script 
rstudioapi::documentSave()

# Current script name:
script_path <- this.path::this.path()

# Load current script
lines <- readLines(script_path, warn = FALSE)

# Find the breakpoint line in current script
startpoint_line <- grep("TEST_STARTPOINT", lines)[1]
breakpoint_line <- grep("TEST_BREAKPOINT", lines)[1]

# Remove all instances of "'PD',"
lines <- gsub("'Bee',", "", lines, fixed = TRUE)

# Replace all instances of another string (example: replace "old" with "new")
lines <- gsub("Manuscript", "Manuscript/dropout", lines, fixed = TRUE)

# Remove everything after breakpoint
truncated_lines <- lines[startpoint_line:breakpoint_line]

source(textConnection(truncated_lines))


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
    save_kable(paste0('Out/Manuscript/tables/alpha_', idx,'_GTDB10.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'KB45_GTDB', 'SM_gtdb-rs220-rep') %>% 
    save_kable(paste0('Out/Manuscript/tables/alpha_', idx,'_GTDB45.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'KB90_GTDB', 'SM_gtdb-rs220-rep') %>% 
    save_kable(paste0('Out/Manuscript/tables/alpha_', idx,'_GTDB90.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'KB10', 'SM_RefSeq_20250528') %>% 
    save_kable(paste0('Out/Manuscript/tables/alpha_', idx,'_RefSeq10.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'KB90', 'SM_RefSeq_20250528') %>% 
    save_kable(paste0('Out/Manuscript/tables/alpha_', idx,'_RefSeq90.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'KB45', 'SM_RefSeq_20250528') %>% 
    save_kable(paste0('Out/Manuscript/tables/alpha_', idx,'_RefSeq45.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'MPA_db2025', 'MOTUS4') %>% 
    save_kable(paste0('Out/Manuscript/tables/alpha_', idx,'_Markers.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'SM_gtdb-rs220-rep', 'SM_RefSeq_20250528') %>% 
    save_kable(paste0('Out/Manuscript/tables/alpha_', idx,'_SM.html'))
  
  quantify_div_variation(
    alpha_div, idx = idx,
    'KB45_GTDB', 'KB45') %>% 
    save_kable(paste0('Out/Manuscript/tables/alpha_', idx,'_KB.html'))
}

# Mean Dataset alphadiv range across methods
# Valid? ??
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



