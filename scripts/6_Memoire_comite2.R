library(pacman)
p_load( magrittr, tidyverse, purrr, kableExtra, phyloseq, patchwork, 
        rstatix, parallel, reshape2, vegan, RColorBrewer)
ps_rare.ls <- read_rds('Out/ps_rare.ls.rds')
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R'))
source("scripts/myFunctions.R")
theme_set(theme_light())

# Work from the species-level table
taxRanks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

################################
# Taxonomic assignment table ####
##################################
# Using tax table only, grouped at each taxrank

tax_assignment <- imap(ps_rare.ls, function(ps_dataset.ls, dataset){
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
}) %>% bind_rows() %>% 
  mutate(Rank = factor(Rank, levels = taxRanks),
         Database = factor(Database, names(tool_colours))) %>% 
  left_join(CCE_metadata, by = 'Database')

# Prepare plot data
these_databases <- c('SM_genbank-2022.03', 'SM_gtdb-rs220-rep', 
                     'MPA_db2022','MPA_db2023', 'MOTUS',
                     'KB10', 'KB45','KB90', 'KB10_GTDB', 'KB45_GTDB','KB90_GTDB')
these_datasets <- c('Moss', 'NAFLD', 'P19_Gut', 'P19_Saliva', 'PD', 'Olive', 'Bee', 'AD_Skin')

# SUBSET
tax_assignment.pdat <- tax_assignment %>% 
  filter(Database %in% these_databases &
           Dataset %in% these_datasets &
           Rank == 'Species') %>% 
  mutate(Prop_db = Num_tax / Num_species_in_db)

#  /$$$$$$$  /$$        /$$$$$$  /$$$$$$$$         /$$  
# | $$__  $$| $$       /$$__  $$|__  $$__/       /$$$$  
# | $$  \ $$| $$      | $$  \ $$   | $$         |_  $$  
# | $$$$$$$/| $$      | $$  | $$   | $$           | $$  
# | $$____/ | $$      | $$  | $$   | $$           | $$  
# | $$      | $$      | $$  | $$   | $$           | $$  
# | $$      | $$$$$$$$|  $$$$$$/   | $$          /$$$$$$
# |__/      |________/ \______/    |__/         |______/

# Number of Species per dataset
tax_assignment.pdat %>% 
  filter(! Database %in% c('SM_gtdb-rs214-rep_MAGs', 'SM_gtdb-rs214-full')) %>% 
  ggplot(aes(y = Dataset, x = Num_tax, fill = Database)) +
  geom_col(position = position_dodge2(preserve = 'single'),
           width = 1.1) + # keep bars the same width
  facet_grid(Dataset~., scales = 'free', switch = 'y') +
  scale_fill_manual(values = tool_colours, labels = CCE_names) +
  labs(y = 'Jeu de données', 
       x = "Nombre d'espèces détectées",
       fill = "Méthodologie") +
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

ggsave('Out/comite2/species_count.pdf',
       bg = 'white', width = 2600, height = 1400,
       units = 'px', dpi = 220)


# Extreme KB10 values
tax_assignment.pdat %>% 
  filter(Tool %in% c('KB10', 'KB10_GTDB')) 

# Min max across all
tax_assignment.pdat %>% 
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
  filter(Database %in% these_databases & Dataset %in% these_datasets) %>% 
  pivot_wider(names_from = Rank, values_from = Num_tax, values_fill = list(Num_tax = 0))

# Split the data by Dataset and create a kable table for each
tables <- tax_assignment_wide %>%
  arrange(CCE_approach) %>% 
  split(.$Dataset) %>%  # Split by Dataset
  imap(function(data, dataset_name) {
    data %>%
      select(-Dataset) %>%  # Remove the Dataset column for the table
      kable(caption = paste("Dataset:", dataset_name), align = "c") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
  })

tables
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

alpha_div <-  imap(ps_rare.ls, function(ps_dataset.ls, dataset){
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
}) %>% list_rbind() %>% 
  # Reorder factors 
  mutate(Database = factor(Database, levels = names(tool_colours))) %>% 
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

these_databases <- c('SM_genbank-2022.03', 'SM_gtdb-rs220-rep', 
                     'MPA_db2022','MPA_db2023', 'MOTUS',
                     'KB10', 'KB10_GTDB',
                     'KB45','KB90', 'KB45_GTDB','KB90_GTDB')
these_datasets <- c('Moss', 'NAFLD', 'P19_Gut', 'P19_Saliva', 'PD', 'Olive', 'Bee', 'AD_Skin')

# Function for a single plot, one plot per CCE approach that includes multiple indices
plot_div_by_approach <- function(df, approach) {
  df %>% filter(approach == CCE_approach) %>% 
  ggplot(aes(x = Database, y = Index_value, )) +
    geom_line(aes(group = Sample), alpha = 0.5, linewidth = 0.2) +
    geom_point(aes(colour = Database)) + 
    facet_grid(Index ~ ., scales = 'free') +
    scale_colour_manual(values = tool_colours, labels = CCE_names) +
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
  filter(Index %in% c('Hill_1')) %>% 
  semi_join(these_samples, join_by(Sample, CCE_approach))%>% 
  plot_div_by_approach(approach = 'DNA-to-Marker')

# PLOT hill_1 variation for methods most equivalent in terms of number of species
these_databases <- c('KB45', 'KB45_GTDB', 'SM_gtdb-rs220-rep', 'SM_genbank-2022.03',
                     'MPA_db2022','MPA_db2023', 'MOTUS')
hill1_comparison <- alpha_div %>% 
  filter(Index %in% c('Hill_1')
          & Database %in% these_databases
           & ! Dataset %in% c('Olive', 'Bee', 'Moss')
         ) %>% 
  mutate(Database = factor(Database, levels = these_databases)) 

hill1_comparison %>% 
  ggplot(aes(x = Database, y = Index_value)) +
  geom_line(aes(group = Sample), alpha = 0.5, linewidth = 0.1) +
  facet_wrap(.~CCE_approach, scales = 'free',
             strip.position = 'bottom') +
  geom_boxplot(aes(fill = Database), alpha = 0.8,
               linewidth = 0.2, outliers = FALSE) + 
  scale_fill_manual(values = tool_colours, labels = CCE_names) +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(fill = 'Méthodologie', x = '', y = "Nombre Hill d'ordre 1")
  
ggsave('Out/comite2/hill1_comparison.pdf',
       bg = 'white', width = 2200, height = 1600,
       units = 'px', dpi = 220)

# For PPT 
ggsave('Out/comite2/hill1_comparison_wide.pdf',
       bg = 'white', width = 2400, height = 1400,
       units = 'px', dpi = 220)


# Quantify those variations
quantify_div_variation <- function(df, ds1, ds2, by_dataset = TRUE) {
  message(paste("Changement entre", ds1, "et", ds2))
  df %>% 
    select(Sample, Database, Index_value, Dataset) %>% 
    filter(Database %in% c(ds1,ds2)) %>% 
    pivot_wider(names_from = Database,
                values_from = Index_value) %>% 
    mutate(change = .[[ds1]] - .[[ds2]]) %>%  # Dynamic column name in mutate()
    filter(!is.na(change)) %>% 
    { if (by_dataset) group_by(., Dataset) else . } %>% 
    summarise(mean_change = mean(change),
              sd_change = sd(change),
              cv = sd_change / mean_change)
}

quantify_div_variation(hill1_comparison, 
                       'KB45_GTDB', 'SM_gtdb-rs220-rep', 
                       by_dataset = TRUE)

quantify_div_variation(hill1_comparison, 
                       'KB45_GTDB', 'KB45', 
                       by_dataset = TRUE)

quantify_div_variation(hill1_comparison, 
                       'SM_gtdb-rs220-rep', 'SM_genbank-2022.03',
                       by_dataset = FALSE)

quantify_div_variation(hill1_comparison,
                       'MPA_db2023', 'MOTUS',
                       by_dataset = FALSE)

# Mean Dataset alphadiv range across methods
alpha_div_subset %>% group_by(Database, Dataset, Taxonomy) %>% 
  filter(Index %in% c('Hill_1') & 
           CCE_approach == 'DNA-to-DNA') %>%
  # Mean tool div by dataset
  summarise(mean = mean(Index_value)) %>% 
  # Mean 
  group_by(Dataset, Taxonomy) %>% 
  summarise(min = min(mean), max = max(mean)) %>% 
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
  filter( #! Dataset %in% c('RA_Gut', 'AD_Skin') &
            Index == 'Shannon') %>% 
  compute_sign_specific_changes('SM_genbank-2022.03', 'SM_gtdb-rs214-rep')

alpha_div %>% 
  filter( ! Dataset %in% c('Moss') &
            Index == 'Shannon' ) %>% 
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
              p.adjust.method = NULL) %>% # conservative
  add_significance() %>% 
  mutate(p.signif = case_when(p.signif == '*' ~ 'p < 0.05',
                              p.signif == '**' ~ 'p < 0.01',
                              p.signif == '***' ~ 'p < 0.001',
                              p.signif == 'ns' ~ 'p > 0.05')) %>% 
  select(Dataset, Database, Index, p, p.signif) %>% 
  mutate(p.signif = factor(p.signif, levels = c('p > 0.05', 'p < 0.05', 'p < 0.01', 'p < 0.001')))
  
alpha_div %>% 
  filter(Index == 'Hill_2' 
        # & Database %in% c('KB10','KB10_GTDB','KB45', 'KB45_GTDB', 'KB90', 'KB90_GTDB', 'MOTUS', 'MPA_db2023', 'SM_gtdb-rs220-rep')
        & ! Dataset %in% c('Bee', 'Moss', 'Olive')
        #& CCE_approach == 'DNA-to-DNA'
        & CCE_approach == 'DNA-to-Marker'
        ) %>% 
  left_join(alpha_div_test, by = c('Dataset', 'Database', 'Index')) %>% 
  ggplot(aes(x = Grouping_var, y = Index_value, fill = p.signif)) +
  geom_boxplot(outliers = FALSE, linewidth =0.2) +
  facet_grid(Dataset~Database, scales = 'free') +
  labs(fill = 'Valeur p', x = 'Groupe', y = 'Nombre de Hill (ordre 1)') +
  scale_y_continuous(limits = c(0, NA)) +
  theme(legend.position = c(0.72,0.95),
        legend.title = element_blank(),
        legend.background = element_rect(
          fill = "white",        # White background
          color = "black",       # Black border
          linewidth = 0.5        # Border thickness
        )) +
  guides(fill = guide_legend(nrow = 1)) 

# Comité ppt
ggsave('Out/comite2/hill1_group_test_wide.pdf',
       bg = 'white', width = 2400, height = 1400,
       units = 'px', dpi = 220)

ggsave('Out/comite2/hill1_group_test_wide2.pdf',
       bg = 'white', width = 2000, height = 2400,
       units = 'px', dpi = 220)

###################
# Beta diversity ###
#####################

# using collapsed pairwise matrices (BC and rAitchison):
# Sample_pair, Dataset, idx, idx_value, CCE, Approach

# pcoa.ls <- read_rds('Out/pcoa_noVST.ls.RDS')
pcoa.ls <- list()
for (dist_metric in c('bray', 'robust.aitchison', 'hellinger')) {
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
  filter(Dist == 'bray' &
           Database %in% these_databases &
           Dataset %in% these_datasets) %>% 
  ggplot(aes(y = Distance, x = Dataset, colour = Database))+
  geom_jitter(size = 0.1, alpha = 0.5,
    position = position_jitterdodge(
      jitter.width = 0.4,  # Adjust jitter width
      dodge.width = 0.7    # Adjust dodge width
  )) +
  geom_boxplot(outliers = FALSE, alpha = 0.7) +
  facet_grid(.~CCE_approach, scale = 'free')+
  scale_colour_manual(values = tool_colours, labels = CCE_names) +
  labs(colour = 'Méthodologie',
       y = 'Dissimilarité bray-curtis', x = 'Jeu de données') 

ggsave(plot = p, 'Out/comite2/change_bc.png', bg = 'white', width = 2600, height = 1400, 
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
database_pairs_of_interest <- list(
  `Kraken 0.45: GTDB – NCBI` = c('KB45_GTDB', 'KB45'),
  `GTDB: Sourmash – Kraken 0.45` = c('SM_gtdb-rs220-rep', 'KB45_GTDB'),
 # `NCBI: Sourmash – Kraken 0.45` = c('SM_genbank-2022.03', 'KB45'),
 `mOTUs – MetaPhlAn 2023` = c('MOTUS', 'MPA_db2023'),
 `MetaPhlAn: 2022 – 2023` = c('MPA_db2022', 'MPA_db2023'),
 `Sourmash: rs220 – rs214` = c('SM_gtdb-rs220-rep', 'SM_gtdb-rs214-rep')
)

# Iterate over pairs of interest
pairwise_dist_gap.df <- imap(database_pairs_of_interest, 
                             function(tool_pair, pair_name){
  compute_distance_differences(pairwise_distances, 
                               tool_pair[1], tool_pair[2]) %>% 
    mutate(
      Pair_name = pair_name
    )
}) %>% list_rbind

# Summary
pairwise_dist_gap.df %>% 
  filter(Dist == 'bray' 
         #& Dataset == 'AD_Skin'
         ) %>% 
  group_by(Pair_name) %>% 
  summarise(mean_diff = mean(abs_diff),
            sd_diff = sd(abs_diff),
            median_diff = median(abs_diff),
            min_diff = min(abs_diff),
            max_diff = max(abs_diff)
  )

# PLOT Boxplot with differences between pairs of interest
pairwise_dist_gap.df %>% 
  filter(Dist == 'bray') %>% 
  mutate(Pair_name = factor(Pair_name, levels = names(database_pairs_of_interest))) %>% 
  ggplot(aes(x = Dataset, y = dist_diff, fill = Dataset)) +
  geom_violin(linewidth = 0.2) +
  facet_grid(.~Pair_name, scale = 'free') +
  theme_bw() +
  labs(y = 'Différences des dissimilarités en paires',
       fill = 'Jeu de données') +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x.top = element_text(angle = 0, hjust = 0),
        panel.grid.major.x = element_blank(),
        legend.position = 'bottom') +
  guides(
    fill = guide_legend(nrow = 1))

ggsave('Out/comite2/bray_diff.pdf', bg = 'white', width = 2200, height = 1200, 
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
  filter(Index == 'hellinger' &
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
                     permutations = 999,
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

# Presented as in poster? 
p_value_lines <- function() {
  list(
    geom_vline(aes(xintercept = log10(0.05), linetype = "p = 0.05"), 
               color = "red", linewidth = 0.3),
    geom_vline(aes(xintercept = log10(0.01), linetype = "p = 0.01"), 
               color = "blue", linewidth = 0.3),
    scale_linetype_manual(name = "Valeur p", 
                          values = c("p = 0.05" = "dashed", 
                                     "p = 0.01" = "dashed"))
  )
}

plot_permanova <- function(ds) {
  ggplot(ds, aes(y = R2, x = log10(p))) +
    geom_point(size = 2, colour = 'black', shape = 23, alpha = 0.8,
               aes(fill = Database)) +
    facet_grid(Dataset~., scales = 'free')  +
    p_value_lines() +
    scale_fill_manual(values = tool_colours, labels = CCE_names) +
    expand_limits(y = 0) +  
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = 'log10(valeur p)', y = '', fill = 'Méthodologie')
}

these_databases <- c('KB10','KB10_GTDB','KB45', 'KB45_GTDB', 'KB90', 'KB90_GTDB', 'MOTUS', 'MPA_db2023','MPA_db2022',  'SM_gtdb-rs220-rep', 'SM_genbank-2022.03')
p_load(patchwork)
p1 <- permanova.ds %>% 
  filter(Index == 'bray'
           & Database %in% these_databases
            & ! Dataset %in% c('Bee', 'Moss', 'Olive', 'AD_Skin')
  ) %>% plot_permanova() +
  labs(y = "Proportion de variance expliquée (R^2)")
  

p2 <- permanova.ds %>% 
  filter(Index == 'bray'
         & Database %in% these_databases
         & Dataset %in% c('Bee', 'Moss', 'Olive', 'AD_Skin')
  ) %>% plot_permanova()

p1 + p2 + plot_layout(guides = 'collect')

ggsave('Out/comite2/permanova_bray.pdf', bg = 'white', width = 2200, height = 1200, 
       units = 'px', dpi = 220)

# descriptive Statistics
permanova.ds %>% 
  filter(Index == 'bray' & Database %in% these_databases) %>% 
  group_by(Dataset) %>% 
  summarise(
    minR2 = min(R2),
    maxR2 = max(R2),
    n = n()
  ) %>% 
  mutate(range = maxR2 - minR2,
         ratio = maxR2/minR2)

########################################################
# Taxa discovery rate (rarefaction curves comparison) ###
##########################################################

