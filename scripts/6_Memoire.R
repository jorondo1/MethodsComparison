library(pacman)
p_load( magrittr, tidyverse, purrr, kableExtra, phyloseq, patchwork, beanplot, 
        rstatix, parallel, reshape2, vegan, RColorBrewer, ggridges, htmltools)
ps_rare.ls <- read_rds('Out/ps_rare.ls.rds')
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R'))
source("scripts/myFunctions.R")
theme_set(theme_light())

# Work from the species-level table
taxRanks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

################################
# Taxonomic assignment table ####
##################################

# Number of taxa by taxrank at the dataset level
tax_assign_ds <- imap(ps_rare.ls, function(ps_dataset.ls, dataset){
  imap(ps_dataset.ls, function(ps, database){
    tax_count <- ps %>% tax_table() %>% data.frame() %>% tibble()
    
    # Iterate over taxRanks
    map_dfr(taxRanks, function(rank) {
      
      num_tax <- tax_count %>%
        filter(!is.na(!!sym(rank))) %>%
        pull(!!sym(rank)) %>%
        unique() %>%
        length()
      
      tibble(
        Dataset = dataset,
        Database = database,
        Rank = rank,
        Num_tax = num_tax
      ) 
    })
  }) %>% bind_rows() 
}) %>% bind_rows()  %>% 
  mutate(Rank = factor(Rank, levels = taxRanks),
         Database = factor(Database, names(tooldb_colours))) %>% 
  left_join(CCE_metadata, by = 'Database')

# Alt: Number of species by sample
tax_assign_sam <- imap(ps_rare.ls, function(ps_dataset.ls, dataset){
  imap(ps_dataset.ls, function(ps, database){
    
    # Species per sample
    psflashmelt(ps) %>% 
      filter(Abundance > 0) %>% 
      select(Sample, Species) %>% 
      distinct() %>% 
      group_by(Sample) %>% 
      summarise(Num_tax = n()) %>% 
      mutate(
        Dataset = dataset,
        Database = database
      )
  }) %>% bind_rows() 
}) %>% bind_rows()  %>% 
  mutate(Database = factor(Database, names(tooldb_colours))) %>% 
  left_join(CCE_metadata, by = 'Database')

# Prepare plot data
these_databases <- c('SM_genbank-2022.03', 'SM_gtdb-rs220-rep', 'SM_RefSeq_20250528',
                     'MPA_db2022','MPA_db2023', 'MOTUS',
                     'KB10', 'KB45','KB90', 'KB10_GTDB', 'KB45_GTDB','KB90_GTDB')
these_datasets <- c('Moss', 'NAFLD', 'P19_Gut', 'P19_Saliva', 'PD', 'Olive', 'Bee', 'AD_Skin')

# SUBSET
tax_assign_ds.pdat <- tax_assign_ds %>% 
  filter(Database %in% these_databases &
           Dataset %in% these_datasets &
           Rank == 'Species') %>% 
  mutate(Prop_db = Num_tax / Num_species_in_db)

tax_assign_sam.pdat <- tax_assign_sam %>% 
  filter(Database %in% these_databases &
           Dataset %in% these_datasets) %>% 
  mutate(Prop_db = Num_tax / Num_species_in_db)

#  /$$$$$$$  /$$        /$$$$$$  /$$$$$$$$         /$$  
# | $$__  $$| $$       /$$__  $$|__  $$__/       /$$$$  
# | $$  \ $$| $$      | $$  \ $$   | $$         |_  $$  
# | $$$$$$$/| $$      | $$  | $$   | $$           | $$  
# | $$____/ | $$      | $$  | $$   | $$           | $$  
# | $$      | $$      | $$  | $$   | $$           | $$  
# | $$      | $$$$$$$$|  $$$$$$/   | $$          /$$$$$$
# |__/      |________/ \______/    |__/         |______/

# Number of Species per datase
tax_assign_ds.pdat %>% 
  ggplot(aes(y = Dataset, x = Num_tax, fill = Database)) +
  geom_col(position = position_dodge2(preserve = 'single',
                                      reverse = TRUE),
           width = 1.1) + # keep bars the same width
  facet_grid(Dataset~., scales = 'free', switch = 'y') +
  scale_fill_manual(values = tooldb_colours, labels = CCE_names) +
  labs(y = 'Dataset', 
       x = "Number of detected species",
       fill = "Methodology") +
  theme(
    panel.spacing.y = unit(0.05, "cm"),
    legend.position = c(.85,.32),
    #axis.text.x = element_text(angle = 45,  hjust=1),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    legend.background = element_rect(
      fill = "white",        # White background
      color = "black",       # Black border
      linewidth = 0.5        # Border thickness
    )
  ) 

ggsave('Out/memoire/species_count_ds.pdf',
       bg = 'white', width = 2600, height = 1400,
       units = 'px', dpi = 220)

# Number of Species per sample per dataset
tax_assign_sam.pdat %>% 
  filter(!str_detect(Database, 'KB10')) %>% 
  ggplot(aes(y = "", x = Num_tax, fill = Database)) +
  geom_boxplot(outliers = FALSE, linewidth = 0.3) + 
  facet_wrap(.~Dataset, scales = 'free', ncol = 2) +
  scale_fill_manual(values = tooldb_colours, labels = CCE_names) +
  labs(y = 'Dataset', 
       x = "Number of detected species",
       fill = "Methodology") +
  theme(
    panel.spacing.y = unit(0.05, "cm"),
    legend.position = c(.85,.15),
    #axis.text.x = element_text(angle = 45,  hjust=1),
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    legend.background = element_rect(
      fill = "white",        # White background
      color = "black",       # Black border
      linewidth = 0.5        # Border thickness
    )
  ) 

ggsave('Out/memoire/species_count_sam.pdf',
       bg = 'white', width = 2600, height = 1400,
       units = 'px', dpi = 220)


# Min max across all
tax_assign_ds.pdat %>% 
  filter(Taxonomy != 'GTDB') %>% 
  group_by(Dataset) %>% 
  summarise(
    min_species = min(Num_tax),
    max_species = max(Num_tax),
    n = n()
  ) %>% 
  mutate(ratio = max_species / min_species) 

## TABLES 
# Reshape the data for each Dataset
tax_assignment_wide <- tax_assignment %>%
  filter(Database %in% these_databases 
         & Dataset %in% these_datasets
         & Rank == 'Species') %>% 
  select(Dataset, Database, Num_tax) %>% 
  pivot_wider(names_from = Dataset, values_from = Num_tax, values_fill = list(Num_tax = 0))

# Split the data by Dataset and create a kable table for each
tax_assignment_wide %>%
  kable( align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>% 
  save_kable(paste0('Out/memoire/tables/species_count.html'))

# By taxrank, number of taxa found by tool, shown independently by approach
# and within approach, intersect size
# One table per dataset
# Is there a clear break across taxranks ?

## Alternatively : series of Venn diagrams ?

####################
# Alpha diversity ###
######################

# For now we only use two-groups variables:
grouping_variable <- c(
  AD_Skin = 'Gender',
  Moss = 'Compartment',
  NAFLD = 'Group',
  P19_Gut = 'diarr',
  P19_Saliva = 'diarr',
  RA_Gut = 'Group',
  PD = 'Group',
  Bee = 'Group',
  Olive = 'Group'
)

# This function will recode them dynamically as A and B, simply 
recode_factor_AB <- function(factor_var) {
  # Ensure variable has exactly 2 levels
  if (nlevels(factor_var) != 2) {
    stop("The input factor must have exactly 2 levels.")
  }
  
  # Recode the factor to "A" and "B" based on level numbers
  factor(as.numeric(factor_var), levels = 1:2, labels = c("A", "B"))
}

alpha_div_tmp <-  imap(ps_rare.ls, function(ps_dataset.ls, dataset){
  imap(ps_dataset.ls, function(ps, database){
    #Estimate indices
    richness <- estimate_diversity(ps, 'Richness')
    shannon <- estimate_diversity(ps, 'Shannon')
    simpson <- estimate_diversity(ps, 'Simpson')
    tail <- estimate_diversity(ps, 'Tail')
    h1 <- estimate_Hill(ps, 1)
    h2 <- estimate_Hill(ps, 2)
    
    # Dataframe with grouping variable
    samdat <- samdat_as_tibble(ps) %>% 
      # Recode the grouping variable as a A/B factor
      mutate(Grouping_var = !!sym(grouping_variable[dataset]) %>% 
               as.factor %>% recode_factor_AB) %>% 
      select(Sample, Grouping_var)
    
    tibble(
      Sample = names(richness),
      Dataset = dataset,
      Database = database,
      Richness = richness,
      Shannon = shannon,
      Hill_1 = h1,
      Hill_2 = h2,
      Tail = tail,
      Simpson = simpson
    ) %>% # Pivot longer for ggplot 
      pivot_longer(cols = c('Richness', 'Shannon', 'Tail', 'Simpson', 'Hill_1', 'Hill_2'),
                   names_to = 'Index',
                   values_to = 'Index_value') %>% 
      left_join( # Add grouping variable
        samdat, by = 'Sample'
      )
    
  }) %>% list_rbind()
}) %>% list_rbind()

# Reorder factors 
alpha_div <- alpha_div_tmp %>% 
  mutate(Database = factor(Database, levels = names(tooldb_colours))) %>% 
  left_join(CCE_metadata, by = 'Database')

### 1. TECHNICAL COMPARISON

#  /$$$$$$$  /$$        /$$$$$$  /$$$$$$$$        /$$$$$$ 
# | $$__  $$| $$       /$$__  $$|__  $$__/       /$$__  $$
# | $$  \ $$| $$      | $$  \ $$   | $$         |__/  \ $$
# | $$$$$$$/| $$      | $$  | $$   | $$           /$$$$$$/
# | $$____/ | $$      | $$  | $$   | $$          /$$____/ 
# | $$      | $$      | $$  | $$   | $$         | $$      
# | $$      | $$$$$$$$|  $$$$$$/   | $$         | $$$$$$$$
# |__/      |________/ \______/    |__/         |________/

these_databases <- c('SM_genbank-2022.03', 'SM_gtdb-rs220-rep', 'SM_RefSeq_20250528',
                     'MPA_db2022','MPA_db2023', 'MOTUS',
                     'KB10', 'KB10_GTDB',
                     'KB45','KB90', 'KB45_GTDB','KB90_GTDB')
these_datasets <- c('Moss', 'NAFLD', 'P19_Gut', 'P19_Saliva', 'PD', 'Olive', 'Bee', 'AD_Skin')

# Function for a single plot, one plot per CCE approach that includes multiple indices
plot_div_by_approach <- function(df, approach) {
  df %>% filter(approach == CCE_approach) %>% 
    ggplot(aes(x = Database, y = Index_value, )) +
    geom_line(aes(group = Sample), alpha = 0.5, linewidth = 0.1) +
    geom_point(aes(colour = Database)) + 
    facet_grid(Index ~ ., scales = 'free') +
    scale_colour_manual(values = tooldb_colours, labels = CCE_names) +
    theme(
      axis.text.x = element_blank()
    ) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(colour = 'Community Composition\nEstimation Tool &\nReference Database')
}

# Subset to databases of interest
alpha_div_subset <- alpha_div %>% 
  filter(#Index %in% c('Shannon') & 
    Dataset %in% these_datasets &
      Database %in% these_databases)

# subset to samples having all tools per approach
these_samples <- alpha_div_subset %>% 
  group_by(Dataset, Index, Sample, CCE_approach) %>% 
  summarise(num = n(), .groups = 'drop') %>% 
  #filter((num == 4 & CCE_approach == 'DNA-to-DNA') | 
  #        (num == 2 & CCE_approach == 'DNA-to-Marker')) %>% 
  dplyr::select(Sample, CCE_approach)

alpha_div_subset %>% # filter only for sample/approach combinations defined
  filter(Index %in% c('Shannon')) %>% 
  semi_join(these_samples, join_by(Sample, CCE_approach))%>% 
  plot_div_by_approach(approach = 'DNA-to-DNA')

alpha_div_subset %>% 
  filter(Index %in% c('Shannon')) %>% 
  semi_join(these_samples, join_by(Sample, CCE_approach)) %>% 
  plot_div_by_approach(approach = 'DNA-to-Marker')

# PLOT hill_1 variation for methods most equivalent in terms of number of species
these_databases <- c('KB45', 'KB45_GTDB', 'SM_gtdb-rs220-rep', 'SM_RefSeq_20250528',
                     'MPA_db2023', 'MOTUS')
alpha_div %>% 
  # We need a composite refdb/
  filter(Index %in% c('Hill_1')
         & Database %in% these_databases
  ) %>% 
  mutate(Database = factor(Database, levels = these_databases)) %>% 
  ggplot(aes(x = Database, y = Index_value)) +
  geom_violin(aes(fill = Tool), 
              linewidth = 0.3) + 
  geom_line(aes(group = Sample), alpha = 0.5, linewidth = 0.1) +
  facet_wrap(~Taxonomy, scales = 'free') +
  scale_fill_manual(values = tool_colours) +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(fill = 'Methodology', x = '', y = "Hill number of order 1 (effective number of equally abundant species)")

ggsave('Out/memoire/hill1_comparison.pdf',
       bg = 'white', width = 2200, height = 1600,
       units = 'px', dpi = 220)

# Quantify those variations
quantify_div_variation <- function(df, ds1, ds2, idx, by_dataset = TRUE) {
  message(paste("Changement entre", ds1, "et", ds2))
  df %>% 
    select(Sample, Database, Index_value, Dataset, Index) %>% 
    filter(Database %in% c(ds1,ds2)
           & Index == idx) %>% 
    pivot_wider(names_from = Database,
                values_from = Index_value) %>% 
    mutate(change = .[[ds1]] - .[[ds2]]) %>%  # Dynamic column name in mutate()
    filter(!is.na(change)) %>% 
    { if (by_dataset) group_by(., Dataset) else . } %>% 
    summarise(mean_change = mean(change),
              sd_change = sd(change),
              cv = sd_change / mean_change)
}

alpha_vars.ls <- list()
for (idx in c('Richness', 'Hill_1', 'Hill_2')) {
  
  alpha_vars.ls[[idx]][['GTDB']] <- quantify_div_variation(
    alpha_div, idx = idx,
    'KB45_GTDB', 'SM_gtdb-rs220-rep',
    by_dataset = TRUE) 
  
  alpha_vars.ls[[idx]][['RefSeq']] <- quantify_div_variation(
    alpha_div, idx = idx,
    'KB45', 'SM_RefSeq_20250528', 
    by_dataset = TRUE)
  
  alpha_vars.ls[[idx]][['Marker']] <- quantify_div_variation(
    alpha_div, idx = idx,
    'MPA_db2023', 'MOTUS',
    by_dataset = TRUE)
}

imap(alpha_vars.ls, function(set.ls, idx) {
  imap(set.ls, function(out_table, comparison) {
    kable(out_table, align = "c") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>% 
      save_kable(paste0('Out/memoire/tables/', idx,'_',comparison,'.html'))
  })
})


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

# Show how samples change between tools
compute_sign_specific_changes <- function(df, tool1, tool2){
  df %>% 
    filter(Database %in% c(tool1, tool2)) %>% 
    select(Sample, Dataset, Index_value, Database) %>% 
    pivot_wider(names_from = Database, values_from = Index_value) %>% 
    mutate(diff = !!sym(tool1) - !!sym(tool2), .keep = 'unused',
           direction = sign(diff)) %>% 
    group_by(Dataset, direction) %>% 
    summarise(mean = mean(diff),
              min = min(diff),
              max = max(diff),
              n = n()) %>% 
    select(-direction)
}

alpha_div %>% 
  filter(Index == 'Hill_1') %>% 
  compute_sign_specific_changes('SM_gtdb-rs220-rep', 'SM_RefSeq_20250528')

alpha_div %>% 
  filter(Index == 'Hill_1') %>% 
  compute_sign_specific_changes('SM_RefSeq_20250528', 'KB45')

alpha_div %>% 
  filter(Index == 'Hill_1') %>% 
  compute_sign_specific_changes('KB45_GTDB', 'SM_gtdb-rs220-rep')

alpha_div %>% 
  filter(Index == 'Shannon' ) %>% 
  compute_sign_specific_changes('MPA_db2023', 'MOTUS')
# Compute mean change? 

### 2. HYPOTHESIS COMPARISON
# Test between groups
# Simple boxplots with pvalues

#  /$$$$$$$  /$$        /$$$$$$  /$$$$$$$$        /$$$$$$ 
# | $$__  $$| $$       /$$__  $$|__  $$__/       /$$__  $$
# | $$  \ $$| $$      | $$  \ $$   | $$         |__/  \ $$
# | $$$$$$$/| $$      | $$  | $$   | $$            /$$$$$/
# | $$____/ | $$      | $$  | $$   | $$           |___  $$
# | $$      | $$      | $$  | $$   | $$          /$$  \ $$
# | $$      | $$$$$$$$|  $$$$$$/   | $$         |  $$$$$$/
# |__/      |________/ \______/    |__/          \______/ 

alpha_div_test <- alpha_div %>% 
  group_by(Dataset, Database, Index) %>% 
  wilcox_test(as.formula("Index_value ~ Grouping_var"),
              p.adjust.method = NULL) %>%  # because we want to see what happens when you do only one
  add_significance()

alpha_div_test %<>% # conservative
  mutate(p.signif = case_when(p < 0.001 ~ 'p < 0.001',
                              p < 0.01 ~ 'p < 0.01',
                              p < 0.05 ~ 'p < 0.05',
                              TRUE ~ 'p ≥ 0.05')) %>% 
  select(Dataset, Database, Index, p, p.signif) %>% 
  mutate(p.signif = factor(p.signif, levels = c('p ≥ 0.05', 'p < 0.05', 'p < 0.01', 'p < 0.001'
  )))

these_databases <- c('KB45', 'KB90', 'SM_RefSeq_20250528', 
                     'KB45_GTDB', 'KB90_GTDB', 'SM_gtdb-rs220-rep',
                     'MPA_db2023','MOTUS')

alpha_labels <- alpha_div_test %>% 
  filter(Index == 'Hill_1' 
         & Database %in% these_databases) %>% 
  mutate(
    Database = factor(Database, levels = these_databases),
    label = paste0('p=',
                   ifelse(p < 0.01, 
                          format(p, scientific = TRUE, digits = 2), 
                          round(p, 2))),            
    x = Inf, y = Inf,                     # Anchor to panel edge
    hjust = 1.1, vjust = 1.5 # Adjust position
  )

alpha_div %>% 
  filter(Index == 'Hill_1' 
         & Database %in% these_databases
         #& CCE_approach == 'DNA-to-DNA'
         #& CCE_approach == 'DNA-to-DNA'
  ) %>% 
  left_join(alpha_div_test, by = c('Dataset', 'Database', 'Index')) %>% 
  mutate(Database = factor(Database, levels = these_databases)) %>% 
  ggplot(aes(x = Grouping_var, y = Index_value, fill = p.signif)) +
  geom_violin(linewidth =0.2, draw_quantiles = c(0.50)) +
  geom_text(
    data = alpha_labels,
    aes(x = x, y = y, label = label, 
        hjust = hjust, vjust = vjust),
    size = 2, color = "black"
  ) + # Show the mean too :
  # stat_summary(fun = mean, geom = "point", size = 2, shape = 3) +
  facet_grid(Dataset~Database, scales = 'free',
             # Relabel facets :
             labeller = labeller(Database = as_labeller(
               setNames(CCE_metadata$MethodName, CCE_metadata$Database)
             ))) +
  labs(fill = 'p-value', x = 'Group', y = 'Hill number of order 1 (effective number of equally abundant species)') +
  scale_y_continuous(limits = c(0, NA)) +
  theme(legend.position = c(0.9,0.75),
        legend.title = element_blank(),
        axis.text = element_text(size = 4),
        legend.background = element_rect(
          fill = "white",        # White background
          color = "black",       # Black border
          linewidth = 0.2        # Border thickness
        )) 

ggsave('Out/memoire/hill1_group_test.pdf',
       bg = 'white', width = 2400, height = 1800,
       units = 'px', dpi = 220)

###################
# Beta diversity ###
#####################

# using collapsed pairwise matrices (BC and rAitchison):
# Sample_pair, Dataset, idx, idx_value, CCE, Approach

# pcoa.ls <- read_rds('Out/pcoa_noVST.ls.RDS')
pcoa.ls <- list()
for (dist_metric in c('bray', 'robust.aitchison')) {
  pcoa.ls[[dist_metric]] <- lapply(ps_rare.ls, function(ps.ls) {
    mclapply(ps.ls, mc.cores = 7, function(ps) {
      compute_pcoa(ps, dist = dist_metric)
    })
  })
}

write_rds(pcoa.ls, 'Out/pcoa_noVST.ls.RDS')

### 1. TECHNICAL COMPARISON
# Essentially the same plot as alphadiv

# Function to extract unique distances for each pair as long df
compile_pair_distances <- function(dist.mx) {
  dist_matrix <- as.matrix(dist.mx)
  
  upper_indices <- which(upper.tri(dist_matrix), arr.ind = TRUE)
  
  data.frame(
    Sample1 = rownames(dist_matrix)[upper_indices[, 1]],
    Sample2 = colnames(dist_matrix)[upper_indices[, 2]],
    Distance = dist_matrix[upper_indices]
  )
}

# Iterate over all pcoa
pairwise_distances <- imap(pcoa.ls, function(dist.ls, dist) {
  imap(dist.ls, function(dataset.ls, dataset) {
    imap(dataset.ls, function(database.ls, database) {
      
      compile_pair_distances(database.ls[['dist.mx']]) %>% 
        mutate(Dist = dist, 
               Dataset = dataset,
               Database = database) 
      
    }) %>% list_rbind
  }) %>% list_rbind
}) %>% list_rbind %>% # Add metadata
  left_join(CCE_metadata, 
            by = 'Database') %>% tibble()

#  /$$$$$$$  /$$        /$$$$$$  /$$$$$$$$       /$$   /$$
# | $$__  $$| $$       /$$__  $$|__  $$__/      | $$  | $$
# | $$  \ $$| $$      | $$  \ $$   | $$         | $$  | $$
# | $$$$$$$/| $$      | $$  | $$   | $$         | $$$$$$$$
# | $$____/ | $$      | $$  | $$   | $$         |_____  $$
# | $$      | $$      | $$  | $$   | $$               | $$
# | $$      | $$$$$$$$|  $$$$$$/   | $$               | $$
# |__/      |________/ \______/    |__/               |__/

# Either a simple dot-plot, but we don't keep track of samples:
p <- pairwise_distances %>% 
  #filter(Dataset != 'PD') %>% 
  filter(Dist == 'bray' 
         & Database %in% these_databases
         & Dataset %in% these_datasets
  ) %>% 
  ggplot(aes(y = Distance, x = "", fill = Database))+
  geom_violin(linewidth = 0.2, draw_quantiles = c(0.5))+
  # geom_jitter(size = 0.1, alpha = 0.5,
  #             position = position_jitterdodge(
  #               jitter.width = 0.4,  # Adjust jitter width
  #               dodge.width = 0.7    # Adjust dodge width
  #             )) +
  #geom_boxplot(outliers = FALSE, alpha = 0.7) +
  facet_grid(.~Dataset, scale = 'free')+
  scale_fill_manual(values = tooldb_colours, labels = CCE_names) +
  labs(colour = 'Methodology',
       y = 'Bray-curtis dissimilarities between each pair of samples', x = 'Dataset'); p

ggsave(plot = p, 'Out/memoire/change_bc.png', bg = 'white', width = 2600, height = 1400, 
       units = 'px', dpi = 220)

compute_distance_differences <- function(df, tool1, tool2) {
  df %>%
    filter(Database %in% c(tool1, tool2)) %>%
    #Create sample pair ID with consistence and uniqueness
    mutate(Pair = ifelse(Sample1 < Sample2, 
                         paste(Sample1, Sample2, sep = "_"),
                         paste(Sample2, Sample1, sep = "_"))) %>%
    select(Database, Dist, Dataset, Pair, Distance) %>%
    # Pivot wide to manually compute difference
    pivot_wider(names_from = Database, values_from = Distance) %>%
    filter(complete.cases(.)) %>%
    mutate(abs_diff = abs(.[[tool1]] - .[[tool2]]),
           dist_diff = .[[tool1]] - .[[tool2]]) %>% # instead of !!sym() , thanks deepseek
    select(-all_of(c(tool1, tool2)))
}

# Define list of pairs of interest, with names used for facet_grid
db_pairs_eval <- list(
  `A. Kraken 0.45 :\nGTDB 220 – RefSeq` = c('KB45_GTDB', 'KB45'),
  `B. Sourmash :\nGTDB 220 – RefSeq` = c('SM_gtdb-rs220-rep', 'SM_RefSeq_20250528'),
  `C. Sourmash GTDB 220 –\n Kraken RefSeq` = c('SM_gtdb-rs220-rep', 'KB45'),
  `D. Kraken GTDB 220 –\n Sourmash RefSeq` = c('KB45_GTDB', 'SM_RefSeq_20250528'),
  `E. GTDB 220 :\nSourmash – Kraken 0.45` = c('SM_gtdb-rs220-rep', 'KB45_GTDB'),
  `F. RefSeq :\nSourmash – Kraken 0.45` = c('SM_RefSeq_20250528', 'KB45'),
  `G. DNA-to-Marker tools :\nmOTUs3 – MetaPhlAn 2023` = c('MOTUS', 'MPA_db2023')
)

db_pairs_ctrl <- list(
  `A. Sourmash\nGTDB220 – GTDB214` = c('SM_gtdb-rs220-rep', 'SM_gtdb-rs214-rep'),
  `B. MetaPhlAn versions\n2023 – 2022` = c('MPA_db2023', 'MPA_db2022'),
  `C. GTDB taxonomy\n214 Full – 214 Reps` = c('SM_gtdb-rs214-full','SM_gtdb-rs214-rep'),
  `D. NCBI taxonomy\nGenbank – RefSeq` = c('SM_genbank-2022.03', 'SM_RefSeq_20250528')
)

# Iterate over pairs of interest
pairwise_dist_gap.df <- imap(
  c(db_pairs_eval,db_pairs_ctrl),
  function(tool_pair, pair_name){
    compute_distance_differences(pairwise_distances, 
                                 tool_pair[1], tool_pair[2]) %>% 
      mutate(
        Pair_name = pair_name
      )
  }) %>% list_rbind %>% 
  mutate(
    Pair_name = factor(Pair_name, 
                       levels = names(c(db_pairs_eval,db_pairs_ctrl))))

pw_dist_gap_eval.df <- pairwise_dist_gap.df %>% 
  filter(Dist == 'bray'
         & Dataset != 'Olive'
         & Pair_name %in% names(db_pairs_eval)
  )

pw_dist_gap_ctrl.df <- pairwise_dist_gap.df %>% 
  filter(Dist == 'bray'
         & Dataset != 'Olive'
         & Pair_name %in% names(db_pairs_ctrl)
  )

# Summary
pairwise_dist_summary <- 
  rbind(pw_dist_gap_eval.df, pw_dist_gap_ctrl.df) %>% 
  group_by(Pair_name, Dataset) %>% 
  summarise(#mean_diff = mean(dist_diff),
    #sd_diff = sd(dist_diff),
    median_diff = median(dist_diff),
    mad_diff = mad(dist_diff),
    min_diff = min(dist_diff),
    max_diff = max(dist_diff)
  ) 

pairwise_dist_summary %>% 
  mutate(Pair_name = str_replace(Pair_name, "\n", " / ")) %>% # prevents markdown pipes from being added
  kable(align = "l") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "responsive")) %>% 
  save_kable(paste0('Out/memoire/tables/beta_diffs.html'))

# visualise cv
pairwise_dist_summary %>% 
  filter(!str_detect(Pair_name,"NCBI taxonomy") 
         & Pair_name %in% names(db_pairs_eval)) %>% 
  ggplot(aes(y = Pair_name)) +
  geom_point(size = 3, aes(x = median_diff, fill = Dataset), shape = 23, colour = 'black') +
  geom_point(size = 5, aes(x = mad_diff,colour = Dataset), shape = 4, stroke = 1)

# Are my sets distributed normally?
pw_dist_gap_eval.df %>% 
  group_by(Dataset, Pair_name) %>% 
  slice_sample(n = 5000, replace = FALSE) %>% 
  shapiro_test(dist_diff) %>%
  ggplot(aes(x = p, y = Pair_name, colour = Dataset)) +
  geom_jitter(width = 0, height = 0.2)
# Mostly no!

# PLOT Boxplot with differences between pairs of interest
plot_dist_gap <- function(df){
  
  ggplot(df, aes(x = Dataset, y = dist_diff, fill = Dataset)) +
    geom_hline(aes(yintercept = 0), 
               color = "darkred", linewidth = 0.5, linetype = 'dashed') +
    geom_violin(linewidth = 0.2, draw_quantiles = c(0.5)) +
    facet_grid(.~Pair_name, scale = 'free') +
    theme_bw() +
    labs(y = 'Same-pair differences in dissimilarities',
         fill = 'Dataset') +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text.x.top = element_text(
            angle = 0, hjust = 0, size = 9),
          panel.grid.major.x = element_blank(),
          legend.position = c(0.5, 0.07),
          legend.title.align = 0.5,
          legend.background = element_rect(
            fill = "white",        # White background
            color = "black",       # Black border
            linewidth = 0.3        # Border thickness
          )) +
    guides(fill = guide_legend(nrow = 1)) 
}
plot_dist_gap(pw_dist_gap_eval.df)
ggsave('Out/memoire/beta_diff_bray.pdf', bg = 'white', width = 2600, height = 1200, 
       units = 'px', dpi = 200)

plot_dist_gap(pw_dist_gap_ctrl.df)
ggsave('Out/memoire/beta_diff_bray_ctrls.pdf', bg = 'white', width = 2500, height = 1200, 
       units = 'px', dpi = 200)

### 2. HYPOTHESIS COMPARISONS
# 2.1. PCoA comparison
pcoa.ds <- imap(pcoa.ls, function(pcoa_ds.ls, dist) {
  imap(pcoa_ds.ls, function(pcoa_db.ls, ds) {
    imap(pcoa_db.ls, function(pcoa, db) {
      
      pcoa$metadata %>%  # create generic grouping var: 
        rownames_to_column('Sample') %>% 
        tibble %>% 
        transmute(Sample = Sample,
                  Grouping_var = factor(!!sym(grouping_variable[[ds]]), labels = c("Group 1", "Group 2")),
                  Index = dist, 
                  Database = db,
                  Dataset = ds,
                  PCo1 = PCo1,
                  PCo2 = PCo2)
      
    }) %>% list_rbind
  }) %>% list_rbind
}) %>% list_rbind

pcoa.ds %>% 
  filter(Index == 'bray' &
           # Dataset %in% "PD" &
           Database %in% c(these_databases)
  ) %>% 
  ggplot(aes(x = PCo1, y = PCo2)) +
  geom_point(aes(colour = Grouping_var), 
             size = 0.5) +
  facet_grid(Dataset ~ Database) +
  stat_ellipse(aes(fill = Grouping_var),
               level=0.95 , geom = "polygon", alpha = 0.18) +
  # theme(legend.position = c(0.47,0.7)) +
  labs(fill = 'Grouping', colour = 'Grouping')

# Procruste comparison of pcoas ?
# permanova.ds <- read_rds('Out/permanova.ds.RDS')

# 2.2. perMANOVA
permanova.ds <- imap(pcoa.ls, function(pcoa_ds.ls, dist) {
  imap(pcoa_ds.ls, function(pcoa_db.ls, ds) {
    imap(pcoa_db.ls, function(pcoa, db) { # iterate over distances
      
      samData <- pcoa$metadata %>% 
        mutate(Grouping_var = factor(!!sym(grouping_variable[[ds]]), 
                                     labels = c("Group 1", "Group 2")))
      
      # permanova
      res <- adonis2(formula = pcoa$dist.mx ~ Grouping_var, 
                     permutations = 9999,
                     data = samData,
                     na.action = na.exclude,
                     parallel = 8)
      
      # parse r2 and p for each explanatory variable 
      tibble(
        Dataset = ds,
        Database = db,       
        Index = dist,
        R2 = res$R2[1],   
        p = res$`Pr(>F)`[1]
      )
      
    }) %>% list_rbind
  }) %>% list_rbind
}) %>% list_rbind %>% 
  left_join(CCE_metadata, by = 'Database')

write_rds(permanova.ds, 'Out/permanova.ds.RDS')

# PLOT variance with p-values
p_load(patchwork)
p_value_lines <- function() {
  list(
    geom_vline(aes(xintercept = log10(0.05), linetype = "p = 0.05"), 
               color = "red", linewidth = 0.2),
    geom_vline(aes(xintercept = log10(0.01), linetype = "p = 0.01"), 
               color = "blue", linewidth = 0.2),
    scale_linetype_manual(name = "p-value thresholds", 
                          values = c("p = 0.05" = "dashed", 
                                     "p = 0.01" = "dashed"))
  )
}

plot_permanova <- function(ds) {
  ggplot(ds, aes(y = R2, x = log10(p))) +
    geom_point(shape = 4,  size = 3, stroke = 0.8,
               aes(colour = Database),
               position = position_jitter(seed = 1, width = 0.01, height = 0)) +
    facet_grid(Dataset~., scales = 'free')  +
    p_value_lines() +
    scale_colour_manual(values = tooldb_colours, labels = CCE_names) +
    expand_limits(y = 0) +  
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = expression("log"[10]*"(p-value)"), y = '', colour = 'Methodology')
}

these_databases <- c('KB45', 'KB90', 'SM_RefSeq_20250528', 
                     'KB45_GTDB', 'KB90_GTDB', 'SM_gtdb-rs220-rep',
                     'MPA_db2023','MOTUS')
p1 <- permanova.ds %>% 
  filter(Index == 'bray'
         & Database %in% these_databases
         & ! Dataset %in% c('Bee', 'Moss', 'Olive', 'AD_Skin')
  ) %>% plot_permanova() +
  labs(y = expression("Proportion of inertia explained (R"^2*")"))


p2 <- permanova.ds %>% 
  filter(Index == 'bray'
         & Database %in% these_databases
         & Dataset %in% c('Bee', 'Moss', 'Olive', 'AD_Skin')
  ) %>% plot_permanova()

p1 + p2 + plot_layout(guides = 'collect')

ggsave('Out/memoire/beta_permanova_bray.pdf', bg = 'white', width = 2200, height = 1200, 
       units = 'px', dpi = 200)

# descriptive Statistics
beta_perm.tbl <- permanova.ds %>% 
  filter(Index == 'bray' & Database %in% these_databases) %>% 
  group_by(Dataset) %>% 
  summarise(
    minR2 = min(R2),
    maxR2 = max(R2),
    n = n()
  ) %>% 
  mutate(range = maxR2 - minR2,
         ratio = maxR2/minR2) %>% 
  kable(align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")); beta_perm.tbl

save_kable(beta_perm.tbl, 'Out/memoire/beta_perm_table.html')
