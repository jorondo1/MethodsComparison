library(pacman)
p_load( magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
        tidyverse, kableExtra, gghalves,
        patchwork,
        rstatix, PairedData)

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
  Richness = 'Richness (number of species detected)',
  Shannon = 'Shannon entropy',
  Hill_1 = 'Hill number of order 1 (effective number of equally abundant species)',
  Simpson = 'Simpson index',
  Hill_2 = 'Inverse Simpson index',
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
    dplyr::select(-n)
  
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
    keep_paired_samples(idx = idx) %>% 
    mutate(
      Facet = case_when(
        Taxonomy == 'Tool-specific' ~ paste0('A. D2M tools (n = ', paired_N, ')'),
        Taxonomy == 'GTDB' ~ paste0('B. D2D tools + GTDB (n = ', paired_N, ')'),
        Taxonomy == 'NCBI' ~ paste0('C. D2D tools + RefSeq (n = ', paired_N, ')'))
    )
}) %>% setNames(names(axis_desc))


# PANEL 1 : distribution of index values ----------------------------------

alpha_distr.plots <- imap(alpha_paired, function(dat, idx) {
  
  dat %>% 
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
    scale_fill_manual(values = tool_colours, breaks = c('MetaPhlAn4', 'mOTUs4', 'Kraken2', 'Sourmash')) +
    theme_light() +
    theme(
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = 'bottom',
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.box.spacing = unit(-0.5, "lines"),
      !!!strip_theme,
      strip.text.x.top = element_text(
        angle = 0, hjust = 0, size = 10)
    ) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(fill = 'Tool', x = '', y = axis_desc[idx])
  
})

# PANEL 2 : distribution of differences -----------------------------------
# Intéressant : encore plus variable/prononcé avec Hill_2

alpha_diff.fun <- function(dat) {
  
  out <- dat %>% 
    # Text for legend titles presented differently:
    mutate(
      Facet = case_when(
        Taxonomy == 'GTDB' ~ 'D2D: Kraken vs. Sourmash (GTDB)',
        Taxonomy == 'Tool-specific' ~ 'D2M: MetaPhlAn4 (2025) vs. mOTUs4',
        Taxonomy == 'NCBI' ~ 'D2D: Kraken vs. Sourmash (RefSeq)'),
      # Reorder facets manually by matching Taxonomy names
      Facet = factor(Facet, levels = Facet[match(c('Tool-specific', 'GTDB', 'NCBI'), Taxonomy)]),
      # All transparent except first level:
      Alpha = ifelse(Taxonomy == 'Tool-specific', 1, 0.7)
    ) %>% 
    dplyr::select(Dataset, Database, Facet, Sample, Index_value, Alpha) %>% 
    
    # Scale indices (median/mad) within dataset/database combo:
    group_by(Dataset, Database) %>% 
    mutate(scaled_index = (Index_value-median(Index_value))/mad(Index_value)) %>% 
    ungroup() %>% 
    
    # Pairwise differences in scaled diversity :
    group_by(Dataset, Facet, Sample, Alpha) %>%
    
    # each group will have 2 by design (see these_databases vector)
    # so we can subtract the first one from the last one within group
    # if they are ordered :
    arrange(Dataset, Facet, Sample, Alpha) %>% 
    summarise(differences = first(scaled_index)-last(scaled_index),
              .groups = 'drop') 
  
  return(out)
  
}

these_indices <- c('Richness', 'Shannon', 'Hill_2', 'Tail')

alpha_diff.pdat <- alpha_paired %>% 
  keep_at(these_indices) %>% 
  map(alpha_diff.fun)

names(alpha_diff.pdat) <- c('Richness', 'Shannon', 'Inverse Simpson', 'Tail')

# --- plotting function
distr_alpha_diffs.fun <- function(plot.dat, plot_name){
  # for anchor poslter figure:
  i      <- match(plot_name, names(alpha_diff.pdat))
  
  plot.dat %>% 
    mutate(strip_label = 'D. Distribution of differences in scaled diversity',
           strip_label_ANCHOR = paste0(LETTERS[i + 3], ". ", plot_name)) %>% 
    # --- TRUNCATE X AXIS for lisibility :
    filter(differences>-2.5, differences < 2.5) %>% 
    ## --- !
    ggplot(aes(x = differences, fill = Facet, alpha = Alpha)) +
    geom_density(linewidth = 0.3) +
    scale_fill_manual(values = c("#7FB364","#823D51","#515A82")) +
    scale_alpha_identity() +
    labs(subtitle = paste0(LETTERS[which(names(alpha_diff.pdat) == plot_name)], ": ", plot_name),
         fill = "Methodology pairs comparison",
         x = 'Differences in scaled diversity', y = 'Density') +
    theme(legend.position = 'bottom',
          panel.grid = element_blank(),
          !!!strip_theme,
          axis.text.y = element_blank()) 
}

distr_alpha_diffs.plots <- alpha_diff.pdat %>% 
  imap(distr_alpha_diffs.fun)  # no change to the call

# Combo plot for paper ------------

# Density distribution plots formatting
distr_alpha_diffs_formatted.plots <- map(distr_alpha_diffs.plots, function(diff_plot){
  
  distr_diffs <- diff_plot +
    facet_wrap(~strip_label) +
    theme(
      plot.subtitle = element_blank(),
      legend.position = c(0.2,0.5),
      strip.text.x.top = element_text(
        angle = 0, hjust = 0, size = 10)
    )    
})

# Half-violin plots formatting
alpha_distr_formatted.plots <- map(alpha_distr.plots, function(distr_plot){
  distr_plot + 
    facet_wrap(~Facet, scales = "free_x") +
    theme(
      axis.ticks.x = element_blank(),
      legend.position = c(0.5,0),
      legend.direction = 'horizontal',
    )
})

# Figure 1 : Shannon plots patchwork
alpha_distr_formatted.plots$Shannon / distr_alpha_diffs_formatted.plots$Shannon +
  plot_layout(design = "
              A
              A
              A
              B
              B") &
  theme(legend.text = element_text(size = 9),
        legend.background = element_rect(
          linewidth = 0.1,
          fill = 'white',
          colour = 'black'))

ggsave(paste0('Out/Manuscript/1.1.alpha_diff_Shannon.pdf'),
       bg = 'white', 
       width = 2250, height = 2571, units = 'px', dpi = 300)

ggsave('Out/Manuscript/PLOS/Fig1.tiff',
       bg = 'white', device = 'tiff', compression = 'lzw',
       width = 2250, height = 2571, units = 'px', dpi = 300)


# Figure SUPP 1 : All distributions, 1 plot each
# No facet grid, allow y axis to be free
imap(alpha_distr.plots, function(distr_plot, idx_name){
  
  ggsave(paste0('Out/Manuscript/1.2.alpha_distr_',idx_name,'.pdf'),
         plot=distr_plot,
         bg = 'white', 
         width = 2250, height = 2200, units = 'px', dpi = 280)
  
})

# Figure SUPP 2 : All difference distributions, in 2 columns -----------
all_distr_diffs <- wrap_plots(distr_alpha_diffs.plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 10),
        strip.text = element_blank(),
        axis.title = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.title = element_blank()) 

ggsave(paste0('Out/Manuscript/1.3.alpha_diff_all.pdf'),
       plot = all_distr_diffs,
       bg = 'white', width = 2100, height = 1400,
       units = 'px', dpi = 200)

# Pitman-Morgan test of differences in variance (paired) ------------------

paired_variance_tests <- function(df, name) {
  
  # Pivot to wide: one row per sample, one column per group
  wide <- df %>%
    dplyr::select(Sample, Facet, differences) %>%
    pivot_wider(
      names_from  = Facet,
      values_from = differences
    )
  
  # Find all pairs:
  groups <- df %>% pull(Facet) %>% unique() %>% sort()
  pairs  <- combn(groups, 2, simplify = FALSE)
  
  # Over all pairs:
  map_dfr(pairs, function(pair) {
    g1 <- as.character(pair[1])
    g2 <- as.character(pair[2])
    
    sub <- wide %>% dplyr::select(all_of(c(g1, g2))) %>% drop_na()
    
    test <- Var.test(sub[[g1]], sub[[g2]], paired = TRUE)
    
    tibble(
      group_1 = g1,
      group_2 = g2,
      n_pairs = nrow(sub),
      var_1   = var(sub[[g1]]),
      var_2   = var(sub[[g2]]),
      var_diff = var_2-var_1,
      t_stat  = unname(test$statistic),
      df      = unname(test$parameter),
      p = test$p.value
    )
  }) %>% 
    mutate(
      Index = name, .keep = 'unused'
    )
}

# Cute table
pitman_all_indices <- imap(alpha_diff.pdat, paired_variance_tests) %>% 
  list_rbind() %>% 
  mutate(p_adj = p.adjust(p, method = 'holm')) 

pitman.table <- pitman_all_indices %>%
  dplyr::select(group_1, group_2, Index, var_1, var_2, var_diff, p_adj) %>% 
  mutate(across(c('var_1', 'var_2', 'var_diff'),~ round(.x, 3)),
         p_adj = round(p_adj, 4)) %T>% 
  print() %>% 
  kable() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

save_kable(pitman.table, 'Out/Manuscript/alpha_pitman.html')


# P values: One pair per row, indices as columns
pitman_all_indices %>%
  dplyr::select(group_1, group_2, p_adj, Index) %>% 
  pivot_wider(
    id_cols = c('group_1', 'group_2'),
    names_from = 'Index',
    values_from = p_adj) %>% 
  kable(format = "simple", digits = 3) %T>% 
  print() %>% 
  writeLines("Out/Manuscript/1.2.alpha_diff_pitmanMorganTest.txt")


# Drop-out test PD --------------------------------------------------------

# >>>> TEST_BREAKPOINT # used to subset the script (stops here)
# save current script 

dropout_fun <- function(dataset){
  rstudioapi::documentSave()
  
  # Current script name:
  script_path <- this.path::this.path()
  
  # Load current script
  lines <- readLines(script_path, warn = FALSE)
  
  # Find the breakpoint line in current script
  startpoint_line <- grep("TEST_STARTPOINT", lines)[1]
  breakpoint_line <- grep("TEST_BREAKPOINT", lines)[1]
  
  # Remove all instances of "'PD',"
  dataset_pattern <- paste0("'", dataset, "',")
  lines <- gsub(dataset_pattern, "", lines, fixed = TRUE)
  
  # Replace all instances of another string (example: replace "old" with "new")
  lines <- gsub("Manuscript", paste0("Manuscript/dropout_",dataset), lines, fixed = TRUE)
  
  # Remove everything after breakpoint
  truncated_lines <- lines[startpoint_line:breakpoint_line]
  
  source(textConnection(truncated_lines), local = TRUE)
}

dropout_fun('Bee')
dropout_fun('PD')


# ANCHOR poster -----------------------------------------------------------

# Density distribution plots formatting
distr_alpha_diffs_formatted_ANCHOR.plots <- map(distr_alpha_diffs.plots, function(diff_plot){
  
  distr_diffs <- diff_plot +
    facet_wrap(~strip_label_ANCHOR) +
    theme(
      legend.position = c(0.2,0.5),
      plot.subtitle = element_blank(),
      strip.background = element_rect(fill = 'grey50'),
      strip.text.x.top = element_text(angle = 0, hjust = 0, size = 12)
    )
})


# Wrap all density plots
all_distr_diffs_ANCHOR <- wrap_plots(distr_alpha_diffs_formatted_ANCHOR.plots, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Distribution of differences in scaled-centered diversity indices.") &
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 10),
        axis.title = element_blank(),
        !!!strip_theme_dark,
        
        panel.spacing = unit(0, "lines"),
        legend.title = element_blank()) 


(alpha_distr_formatted.plots$Shannon + theme(!!!strip_theme_dark)) / 
  wrap_elements(full = all_distr_diffs_ANCHOR) +
  plot_layout(
    design = "
              A
              A
              B
              B") &
  theme(
    legend.text = element_text(size = 10),
    legend.background = element_rect(
      linewidth = 0.1,
      fill = 'white',
      colour = 'black'))

ggsave(paste0('Out/Manuscript/1.4.alpha_diff_Shannon_ANCHOR.pdf'),
       bg = 'white', width = 2300, height = 2400,
       units = 'px', dpi = 220)


# SAAAAAAAAAAANDBOX 

# Spearman correlation of ranks -------------------------------------------

spearman.raw <- imap(alpha_paired, function(dat, idx) {
  
  dat %>% 
    dplyr::select(Taxonomy, Database, Sample, Index_value) %>% 
    group_split(Taxonomy) %>%
    map_dfr(function(x) {
      x %>%
        pivot_wider(names_from = Database, values_from = Index_value) %>%
        dplyr::select(-Sample, -Taxonomy) %>% 
        cor_test(vars = everything(), method = "spearman")  %>%
        transmute(
          Taxonomy = unique(x$Taxonomy),
          Index = idx,
          cor = cor,
          p = p)
    })
}) %>% list_rbind() 

spearman.table <- 
  spearman.raw %>% 
  # hill numbers are monotonic transformations, hence redundant
  filter(!Index %in% c("Hill_1", "Hill_2")) %>% 
  mutate(
    Approach = case_when(
      Taxonomy == 'Tool-specific' ~ paste0('MetaPhlAn4 (2025) vs. mOTUs4'),
      Taxonomy == 'GTDB' ~ paste0('Kraken vs. Sourmash (GTDB)'),
      Taxonomy == 'NCBI' ~ paste0('Kraken vs. Sourmash (RefSeq)')),
    .keep = 'unused', .after = Index
  ) %>% 
  pivot_wider(names_from = Index, values_from = cor, id_cols = Approach) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

save_kable(spearman.table, 'Out/Manuscript/alpha_corr.html')

# Sample diversity variations tables --------------------------------------

quantify_div_variation <- function(df, ds1, ds2, idx) {
  dir_change <- df %>% 
    dplyr::select(Sample, Database, Index_value, Dataset, Index) %>% 
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

for (idx in c('Richness', 'Shannon', 'Hill_2', 'Tail')) {
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
# # Valid? ??
# alpha_div %>% group_by(Database, Dataset, Taxonomy) %>% 
#   filter(Index %in% c('Hill_1')
#          & Database %in% these_databases) %>%
#   # Mean tool div by dataset
#   summarise(mean = mean(Index_value), .groups = 'drop') %>% 
#   # Mean 
#   group_by(Dataset, Taxonomy) %>% 
#   summarise(min = min(mean), max = max(mean), .groups = 'drop') %>% 
#   group_by(Taxonomy) %>% 
#   mutate(fold_increase = max / min) %>% 
#   summarise(mean_fold = mean(fold_increase),
#             sd_fold = sd(fold_increase))
# 
# 
# 
