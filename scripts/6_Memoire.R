library(pacman)
p_load( magrittr, tidyverse, purrr, kableExtra, phyloseq, patchwork, rstatix, parallel, reshape2)
ps_rare.ls <- read_rds('Out/ps_rare.ls.rds')
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R'))
source("scripts/myFunctions.R")
theme_set(theme_light())

# Work from the species-level table
ps_rare.ls <- ps_rare.ls$Species
taxRanks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

################################
# Taxonomic assignment table ####
##################################
# Using tax table only, grouped at each taxrank

tax_assignment <- imap(ps_rare.ls, function(ps_dataset.ls, dataset){
  imap(ps_dataset.ls, function(ps, CCE_tool){
    tax_count <- ps %>% tax_table %>% data.frame %>% tibble
    
    # Iterate over taxRanks
    map_dfr(taxRanks, function(rank) {
      num_tax <- tax_count %>%
        filter(!is.na(!!sym(rank))) %>%
        pull(!!sym(rank)) %>%
        unique() %>%
        length()
      
      tibble(
        Dataset = dataset,
        CCE_tool = CCE_tool,
        Rank = rank,
        Num_tax = num_tax
      ) %>% left_join(select(CCE_metadata, database, CCE_approach), 
                      by = join_by(CCE_tool == database)) # Add approach
    })
  }) %>% bind_rows() 
}) %>% bind_rows() %>% 
  mutate(Rank = factor(Rank, levels = taxRanks),
         CCE_tool = factor(CCE_tool, names(tool_colours)))

# Prepare plot data
these_databases <- c('SM_genbank-2022.03', 'SM_gtdb-rs214-full', 'MPA_db2023',  'KB45','KB90', 'MOTUS')
these_datasets <- c('AD_Skin', 'Moss', 'NAFLD', 'P19_Gut', 'P19_Saliva')

tax_assignment.pdat <- tax_assignment %>% 
  filter(CCE_tool %in% these_databases &
           Dataset %in% these_datasets &
           Rank %in% c('Species', 'Genus', 'Family', 'Order', 'Phylum'))

#  /$$$$$$$  /$$        /$$$$$$  /$$$$$$$$         /$$  
# | $$__  $$| $$       /$$__  $$|__  $$__/       /$$$$  
# | $$  \ $$| $$      | $$  \ $$   | $$         |_  $$  
# | $$$$$$$/| $$      | $$  | $$   | $$           | $$  
# | $$____/ | $$      | $$  | $$   | $$           | $$  
# | $$      | $$      | $$  | $$   | $$           | $$  
# | $$      | $$$$$$$$|  $$$$$$/   | $$          /$$$$$$
# |__/      |________/ \______/    |__/         |______/

tax_assignment.pdat %>% 
  ggplot(aes(x = Dataset, y = Num_tax, fill = CCE_tool)) +
  geom_col(position = position_dodge2(preserve = 'single')) + # keep bars the same width
  facet_grid(Rank~., scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  scale_fill_manual(values = tool_colours, labels = CCE_names) 
  
# Stats
tax_assignment.pdat %>% 
  filter(Rank == 'Species') %>% 
  group_by(Dataset) %>% 
  summarise(min_count = min(Num_tax),
            max_count = max(Num_tax))

# Reshape the data for each Dataset
tax_assignment_wide <- tax_assignment %>%
  filter(CCE_tool %in% these_databases & Dataset %in% these_datasets) %>% 
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
  RA_Gut = 'Group'
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
  imap(ps_dataset.ls, function(ps, CCE_tool){
    
    #Estimate indices
    richness <- estimate_diversity(ps, 'Richness')
    shannon <- estimate_diversity(ps, 'Shannon')
    simpson <- estimate_diversity(ps, 'Simpson')
    tail <- estimate_diversity(ps, 'Tail')
    
    # Dataframe with grouping variable
    samdat <- samdat_as_tibble(ps) %>% 
      # Recode the grouping variable as a A/B factor
      mutate(Grouping_var = !!sym(grouping_variable[dataset]) %>% 
               as.factor %>% recode_factor_AB) %>% 
      select(Sample, Grouping_var)
  
    # Pull CCE approach
    approach <- CCE_metadata %>% # defined in myFunctions.R
      filter(database == CCE_tool) %>% 
      pull(CCE_approach)
    
    # Compile data
    tibble(
      Sample = names(richness),
      Dataset = dataset,
      CCE_tool = CCE_tool,
      CCE_approach = approach,
      Richness = richness,
      Shannon = shannon,
      Tail = tail,
      Simpson = simpson
    ) %>% # Pivot longer for ggplot 
      pivot_longer(cols = c('Richness', 'Shannon', 'Tail', 'Simpson'),
                   names_to = 'Index',
                   values_to = 'Index_value') %>% 
      left_join( # Add grouping variable
        samdat, by = 'Sample'
      )
      
  }) %>% list_rbind()
}) %>% list_rbind() %>% 
  # Reorder factors 
  mutate(CCE_tool = factor(CCE_tool, levels = names(tool_colours)))

### 1. TECHNICAL COMPARISON

#  /$$$$$$$  /$$        /$$$$$$  /$$$$$$$$        /$$$$$$ 
# | $$__  $$| $$       /$$__  $$|__  $$__/       /$$__  $$
# | $$  \ $$| $$      | $$  \ $$   | $$         |__/  \ $$
# | $$$$$$$/| $$      | $$  | $$   | $$           /$$$$$$/
# | $$____/ | $$      | $$  | $$   | $$          /$$____/ 
# | $$      | $$      | $$  | $$   | $$         | $$      
# | $$      | $$$$$$$$|  $$$$$$/   | $$         | $$$$$$$$
# |__/      |________/ \______/    |__/         |________/

# Function for a single plot, one plot per CCE approach that includes multiple indices
plot_div_by_approach <- function(df, approach) {
  df %>% filter(approach == CCE_approach) %>% 
  ggplot(aes(x = CCE_tool, y = Index_value, )) +
    geom_line(aes(group = Sample), alpha = 0.5, linewidth = 0.2) +
    geom_point(aes(colour = CCE_tool)) + 
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
  filter(Index %in% c('Richness', 'Shannon') & 
           Dataset %in% these_datasets &
           CCE_tool %in% these_databases)

# subset to samples having all tools per approach
these_samples <- alpha_div_subset %>% 
  group_by(Dataset, Index, Sample, CCE_approach) %>% 
  summarise(num = n(), .groups = 'drop') %>% 
  filter((num == 4 & CCE_approach == 'DNA-to-DNA') | 
           (num == 2 & CCE_approach == 'DNA-to-Marker')) %>% 
  dplyr::select(Sample, CCE_approach)

alpha_div_subset %>% # filter only for sample/approach combinations defined
  semi_join(these_samples, join_by(Sample, CCE_approach))%>% 
  plot_div_by_approach(approach = 'DNA-to-DNA')

alpha_div_subset %>% 
  semi_join(these_samples, join_by(Sample, CCE_approach))%>% 
  plot_div_by_approach(approach = 'DNA-to-Marker')

# Or One plot per index showing tools from both approaches 
plot_div_by_index <- function(df, index) {
  df %>%   filter(Index == index) %>% 
  ggplot(aes(x = CCE_tool, y = Index_value, )) +
    geom_line(aes(group = Sample), alpha = 0.5, linewidth = 0.2) +
    geom_point(aes(colour = Dataset)) + 
    facet_grid(. ~ CCE_approach, scales = 'free') +
    scale_x_discrete(labels = CCE_names) +
    theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = 'CCE Tool', y = index,
         colour = 'Community Composition\nEstimation Tool &\nReference Database')
}

alpha_div_subset2 <- alpha_div %>% 
  semi_join(these_samples, join_by(Sample, CCE_approach))%>% 
  filter(Dataset %in% these_datasets &
           CCE_tool %in% these_databases) %>% 
  mutate(CCE_tool = factor(CCE_tool, levels = c('SM_genbank-2022.03', 'SM_gtdb-rs214-full', 'KB45','KB90', 'MPA_db2023', 'MOTUS'))) 

plot_div_by_index(alpha_div_subset2, 'Richness')
ggsave('Out/memoire/change_richness.png', bg = 'white', width = 2000, height = 2000, 
       units = 'px', dpi = 220)

plot_div_by_index(alpha_div_subset2, 'Shannon')
ggsave('Out/memoire/change_shannon.png', bg = 'white', width = 2000, height = 2000, 
       units = 'px', dpi = 220)

# Show how samples change between tools
compute_sign_specific_changes <- function(df, tool1, tool2){
  df %>% 
    filter(CCE_tool %in% c(tool1, tool2)) %>% 
    select(Sample, Dataset, Index_value, CCE_tool) %>% 
    pivot_wider(names_from = CCE_tool, values_from = Index_value) %>% 
    mutate(diff = !!sym(tool1) - !!sym(tool2), .keep = 'unused',
           direction = sign(diff)) %>% 
    group_by(Dataset, direction) %>% 
    summarise(mean = mean(diff),
              min = min(diff),
              max = max(diff),
              n = n())
}

alpha_div %>% 
  filter( ! Dataset %in% 'RA_Gut' &
            Index == 'Richness') %>% 
  compute_sign_specific_changes('SM_gtdb-rs214-full', 'KB45')
  
alpha_div %>% 
  filter( ! Dataset %in% 'Moss' &
            Index == 'Shannon') %>% 
  compute_sign_specific_changes('KB90', 'KB45')

alpha_div %>% 
  filter( ! Dataset %in% c('RA_Gut', 'AD_Skin') &
            Index == 'Richness') %>% 
  compute_sign_specific_changes('SM_genbank-2022.03', 'SM_gtdb-rs214-full')

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
  group_by(Dataset, CCE_tool, Index) %>% 
  wilcox_test(as.formula("Index_value ~ Grouping_var"),
              p.adjust.method = NULL) %>% # conservative
  add_significance() %>% 
  select(Dataset, CCE_tool, Index, p, p.signif)
  
alpha_div %>% 
  filter(Index == 'Shannon' 
         & CCE_tool %in% these_databases
         ) %>% 
  left_join(alpha_div_test, by = c('Dataset', 'CCE_tool', 'Index')) %>% 
  ggplot(aes(x = Grouping_var, y = Index_value, colour = p.signif)) +
  geom_boxplot(outliers = FALSE) +
  facet_grid(Dataset~CCE_tool, scales = 'free') +
  scale_y_continuous(limits = c(0, NA)) 

###################
# Beta diversity ###
#####################

# using collapsed pairwise matrices (BC and rAitchison):
# Sample_pair, Dataset, idx, idx_value, CCE, Approach

pcoa.ls <- list()
for (dist_metric in c('bray', 'robust.aitchison')) {
  pcoa.ls[[dist_metric]] <- lapply(ps_rare.ls, function(ps.ls) {
    mclapply(ps.ls, mc.cores = 7, function(ps) {
      compute_pcoa(ps, dist = dist_metric)
    })
  })
}

### 1. TECHNICAL COMPARISON
# Essentially the same plot as alphadiv

# Function to extract unique distances for each pair as long df
compile_pair_distances <- function(dist.mx) {
  
  dist_matrix <- as.matrix(dist.mx)
  dist_matrix[lower.tri(dist_matrix, diag = TRUE)] <- NA #lower triangle
  
  as.data.frame(dist_matrix) %>% # 
    rownames_to_column(var = "Sample1") %>%
    pivot_longer(-Sample1, names_to = "Sample2", values_to = "Distance", values_drop_na = TRUE)
}

pairwise_distances <- imap(pcoa.ls, function(dist.ls, dist) {
  imap(dist.ls, function(dataset.ls, dataset) {
    imap(dataset.ls, function(database.ls, database) {
      
      compile_pair_distances(database.ls[['dist.mx']]) %>% 
        mutate(Dist = dist, 
               Dataset = dataset,
               CCE_tool = database)
    }) %>% list_rbind
  }) %>% list_rbind
}) %>% list_rbind %>% 
  left_join(select(CCE_metadata, database, CCE_approach), 
            by = join_by(CCE_tool == database))

#  /$$$$$$$  /$$        /$$$$$$  /$$$$$$$$       /$$   /$$
# | $$__  $$| $$       /$$__  $$|__  $$__/      | $$  | $$
# | $$  \ $$| $$      | $$  \ $$   | $$         | $$  | $$
# | $$$$$$$/| $$      | $$  | $$   | $$         | $$$$$$$$
# | $$____/ | $$      | $$  | $$   | $$         |_____  $$
# | $$      | $$      | $$  | $$   | $$               | $$
# | $$      | $$$$$$$$|  $$$$$$/   | $$               | $$
# |__/      |________/ \______/    |__/               |__/

# Either a dotplot
pairwise_distances %>% 
  filter(#Dist == 'bray' &
           CCE_tool %in% these_databases &
           Dataset %in% these_datasets) %>% 
  ggplot(aes(y = Distance, x = Dataset, colour = CCE_tool))+
  geom_jitter(size = 0.1,alpha = 0.5,
    position = position_jitterdodge(
      jitter.width = 0.4,  # Adjust jitter width
      dodge.width = 0.7    # Adjust dodge width
  )) +
  geom_boxplot(outliers = FALSE, alpha = 0.7) +
  facet_grid(Dist~CCE_approach, scale = 'free')+
  scale_colour_manual(values = tool_colours, labels = CCE_names)
  
ggsave('Out/memoire/change_beta.png', bg = 'white', width = 2000, height = 2000, 
       units = 'px', dpi = 220)

pairwise_dist_split <- pairwise_distances %>% 
  filter(CCE_tool %in% c('KB45', 'MOTUS', 'SM_genbank-2022.03', 'MPA_db2023')) %>% 
  mutate(Pair = paste0(Sample1, '_',Sample2), .keep = 'unused') %>% 
  group_by(CCE_approach, Dataset, Dist, Pair) %>% 
  select(-CCE_tool) %>% 
  group_split %>% # subset by group (creates a list of tibbles)
  .[sapply(., nrow) == 2] # drop groups that only have one tool

pivot_group <- function(group) {
  group %>% 
    mutate(CCE_ID = paste0("CCE_", row_number())) %>% 
    pivot_wider(names_from = CCE_ID, values_from = Distance) %>% 
    mutate(dist_diff = abs(CCE_1 - CCE_2), .keep = 'unused') 
}

pairwise_dist_gap <- pairwise_dist_split %>% 
  mclapply(pivot_group, mc.cores = 8) %>% 
  bind_rows()

pairwise_dist_gap %>% 
  ggplot(aes(x = Dataset, y = dist_diff, fill = Dataset)) +
  geom_boxplot() +
  # geom_jitter(size = 0.3,
  #             position = position_jitterdodge(
  #               jitter.width = 0.4,  # Adjust jitter width
  #               dodge.width = 0.7    # Adjust dodge width
  #             )) +
  facet_grid(Dist~CCE_approach, scale = 'free') +
  theme_light() +
  labs(y = 'Differences in pairwise distances between tools')

ggsave('test_dist.png', bg = 'white', width = 2000, height = 2000, 
       units = 'px', dpi = 220)


# Find every testable pair of tools from the data 
grouped_pairs <- pairwise_distances %>% 
  select(CCE_tool, CCE_approach, Dataset) %>% 
  group_by(CCE_approach, Dataset) %>% 
  distinct() %>% 
  summarise(Tool = list(expand.grid(CCE_tool, CCE_tool)), .groups = 'drop') %>% 
  unnest(Tool) %>% 
  filter(Var1 != Var2) %>% # remove same tools
  rowwise() %>% 
  mutate(Tool = list(sort(c(Var1, Var2)))) %>%  # Sort pairs to ensure consistency
  distinct(CCE_approach, Dataset, Tool) %>%   # Remove duplicate pairs
  ungroup() %>% 
  unnest_wider(Tool, names_sep = '_') %>% # normal tibble
  mutate(across(everything(), as.character)) 

# Function runs paired wilcox for a given pair of tools 
pairwise_test_wilcox_custom <- function(
    dist_data, ds, dist, tool1, tool2
    ){
  pairwise_wide <- dist_data %>% 
    filter(Dist== dist & # Filter the dataset to subset of interest
             Dataset == ds &
             CCE_tool %in% c(tool1, tool2)) %>% 
    rowwise() %>% # then create sample pair names 
    mutate(Pair = paste0(min(Sample1, Sample2), '_', max(Sample1, Sample2))) %>% 
    pivot_wider(names_from = CCE_tool, values_from = Distance) %>% # one row per pair
    drop_na() # keep only distances that were eval by 2 tools
  
    # Wilcox test
    w.test <- wilcox.test(pull(pairwise_wide,tool1), 
                pull(pairwise_wide,tool2), 
                paired = TRUE)
    
    # Rank-biserial correlation
    rrs <- rank_biserial(pull(pairwise_wide,tool1), 
                         pull(pairwise_wide,tool2), 
                         paired = TRUE)
    
    # Output a tibble
    tibble(Tool1 = tool1,
           Tool2 = tool2,
           Dist = dist, 
           Dataset = ds,
           V = w.test$statistic,
           p = w.test$p.value,
           rrs = rrs$r_rank_biserial)
}

# Iterate the function over every possible pairs
exclude_these_tools <- c('KB51', 'MPA_db2022')

results <- grouped_pairs %>% 
  filter(!Tool_1 %in% exclude_these_tools &
           !Tool_2 %in% exclude_these_tools) %>% 
  # Iteration: 
  pmap(function(CCE_approach, Dataset, Tool_1, Tool_2) {
    pairwise_test_wilcox_custom(
      dist_data = pairwise_distances,
      ds = Dataset, dist = 'bray',
      tool1 = Tool_1, tool2 = Tool_2)
}) %>% list_rbind %>% 
  mutate(p.adj = round(p.adjust(p), 4)) # adjust pvals

### 2. HYPOTHESIS COMPARISONS
# 2.1. PCoA comparison
# Procruste comparison of pcoas ?

# 2.2. perMANOVA 
# Presented as in poster? 

########################################################
# Taxa discovery rate (rarefaction curves comparison) ###
##########################################################

