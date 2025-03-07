library(pacman)
p_load( magrittr, tidyverse, purrr, kableExtra)
ps_rare.ls <- read_rds('Out/ps_rare.ls.rds')
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R'))
source("scripts/myFunctions.R")

# Work from the species-level table
ps_rare.ls <- ps_rare.ls$Species
taxRanks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
which_databases <- c('SM_genbank-2022.03', 'MPA_db2023', 'KB51', 'MOTUS')
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
  mutate(Rank = factor(Rank, levels = taxRanks))

tax_assignment %<>% 
  filter(CCE_tool %in% which_databases &
           Rank %in% c('Species', 'Genus', 'Family')) 

tax_assignment %>% 
  #filter(CCE_tool %in% c('SM_genbank-2022.03', 'MPA_db2023', 'KB20', 'MOTUS')) %>% 
  ggplot(aes(x = Dataset, y = Num_tax, fill = CCE_tool)) +
  geom_col(position = position_dodge()) +
  facet_grid(Rank~CCE_approach, scales = 'free_y') +
  theme_light() +
  scale_fill_manual(values = tool_colours)

# Reshape the data for each Dataset
tax_assignment_wide <- tax_assignment %>%
  filter(CCE_tool %in% which_databases) %>% 
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
  P19_Saliva = 'diarr'
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

# Function for a single plot, one plot per CCE approach that includes multiple indices
plot_div_by_approach <- function(df, approach) {
  df %>% filter(approach == CCE_approach) %>% 
  ggplot(aes(x = CCE_tool, y = Index_value, )) +
    geom_line(aes(group = Sample), alpha = 0.5, linewidth = 0.2) +
    geom_point(aes(colour = CCE_tool)) + 
    facet_grid(Index ~ Dataset, scales = 'free') +
    scale_colour_manual(values = tool_colours, labels = CCE_names) +
    theme_light() +
    theme(
      axis.text.x = element_blank()
    ) +
    labs(colour = 'Community Composition\nEstimation Tool &\nReference Database')
}

alpha_div_subset <- alpha_div %>% 
  filter(Index %in% c('Richness', 'Shannon') & 
           CCE_tool %in% which_databases)

p1 <- alpha_div_subset %>% 
  plot_div_by_approach(approach = 'DNA-to-DNA')
p2 <- alpha_div_subset %>% 
  plot_div_by_approach(approach = 'DNA-to-Marker')

# Or One plot per index showing tools from both approaches 

plot_div_by_index <- function(df) {
  ggplot(df, aes(x = CCE_tool, y = Index_value, )) +
    geom_line(aes(group = Sample), alpha = 0.5, linewidth = 0.2) +
    geom_point(aes(colour = CCE_tool)) + 
    facet_grid(Dataset ~ CCE_approach, scales = 'free') +
    scale_colour_manual(values = tool_colours, labels = CCE_names) +
    theme_light() +
    theme(
      axis.text.x = element_blank()
    ) +
    labs(colour = 'Community Composition\nEstimation Tool &\nReference Database')
}

alpha_div %>% 
  filter(Index == 'Simpson' & 
           CCE_tool %in% which_databases) %>% 
  plot_div_by_index()

# Compute mean change? 

### 2. HYPOTHESIS COMPARISON
# Test between groups
# Simple boxplots with pvalues

alpha_div_test <- alpha_div %>% 
  group_by(Dataset, CCE_tool, Index) %>% 
  wilcox_test(as.formula("Index_value ~ Grouping_var"),
              p.adjust.method = NULL) %>% # conservative
  add_significance() %>% 
  select(Dataset, CCE_tool, Index, p, p.signif)
  

alpha_div %>% 
  filter(Index == 'Shannon' 
         & CCE_tool %in% which_databases
         ) %>% 
  left_join(alpha_div_test, by = c('Dataset', 'CCE_tool', 'Index')) %>% 
  ggplot(aes(x = Grouping_var, y = Index_value, fill = p.signif)) +
  geom_boxplot(outliers = FALSE) +
  facet_grid(Dataset~CCE_tool, scales = 'free') +
  theme_light()

###################
# Beta diversity ###
#####################

# using collapsed pairwise matrices (BC and rAitchison):
# Sample_pair, Dataset, idx, idx_value, CCE, Approach

### 1. TECHNICAL COMPARISON
# Essentially the same plot as alphadiv

### 2. HYPOTHESIS COMPARISONS
# 2.1. PCoA comparison
# Procruste comparison of pcoas ?

# 2.2. perMANOVA 
# Presented as in poster? 



