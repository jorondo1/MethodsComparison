library(pacman)
p_load( magrittr, tidyverse, purrr, kableExtra)
ps_rare.ls <- read_rds('Out/ps_rare.ls.rds')
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R'))
source("scripts/myFunctions.R")

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
  mutate(Rank = factor(Rank, levels = taxRanks))

tax_assignment %>% 
  filter(!CCE_tool %in% c('KB05','SM_gtdb-rs214-rep_MAGs', 'MPA_db2019')) %>% 
  #filter(CCE_tool %in% c('SM_genbank-2022.03', 'MPA_db2023', 'KB20', 'MOTUS')) %>% 
  ggplot(aes(x = Dataset, y = Num_tax, fill = CCE_tool)) +
  geom_col(position = position_dodge()) +
  facet_grid(Rank~CCE_approach, scales = 'free') +
  theme_light() +
  scale_fill_manual(values = tool_colours)

# Reshape the data for each Dataset
tax_assignment_wide <- tax_assignment %>%
  filter(CCE_tool %in% c('SM_genbank-2022.03', 'MPA_db2023', 'KB20', 'MOTUS')) %>% 
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


# By taxrank, number of taxa found by tool, shown independently by approach
# and within approach, intersect size
# One table per dataset
# Is there a clear break across taxranks ?

## Alternatively : series of Venn diagrams ?

####################
# Alpha diversity ###
######################
# Sample_id, Dataset, idx, idx_value, CCE, Approach
# Import from previous work or re-generate?

# Tests at lowest (mostly) coherent taxonomic rank, than one step lower
# Or straight up at species level ? because that's one reason to do shotgun

### 1. TECHNICAL COMPARISON
# One plot per index (Richness and Shannon, to begin with)
# facet dataset ~ approach; x = tool and y = idx_value
# dots + grouped lines
# Compute mean change? 

### 2. HYPOTHESIS COMPARISON
# Test between groups
# Simple boxplots with pvalues

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



