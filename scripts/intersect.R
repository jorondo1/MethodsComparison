# intersect in taxonomy
# pairwise Jaccard index on detected species -> similarity matrix
# First test at genus level
library(pacman)
p_load(magrittr, tidyverse, phyloseq)

# 1. Extract taxonomic tables + convert to P/A

ps_genus.ls <- read_rds("Out/ps_rare_genus.ls.rds") 
ps_family.ls <- read_rds("Out/ps_rare_family.ls.rds") 

# look at names of taxa part of the exclusion space
extract_taxnames <- function(ps, ds, db1, db2, taxRank) {
  
   extract_tax <- function(db) {
     ps[[ds]][[db]] %>% psmelt %>% 
       group_by({{ taxRank }}) %>%
       summarise(Abundance = sum(Abundance)) %>% 
       pull({{ taxRank }}) %>% unique
   }
  tax1 <- extract_tax(db1)
  tax2 <- extract_tax(db2)
  intersect_set <- intersect(tax1,tax2)
  message(paste(length(intersect_set), 'found in common'))
  
  union_set <- union(tax1, tax2)
  message(paste(length(union_set), 'found in total'))
  message(paste(length(intersect_set)/length(union_set), 'intersect'))
  
  setdiff(union_set, intersect_set) %>% sort
}

extract_taxnames(ps_family.ls,
                 'Saliva', 
                 'KB51', 'SM_gtdb_rs214_full', 
                 Family)

########################
# Compute intersect ###
######################

# 

extract_sample_PA <- function(ps_list, taxRank) {
  map(names(ps_list), function(ds){ # first list level: apply to all dataset
    ds_sublist <- ps_list[[ds]] # extract sublist (ps object)
    map(names(ds_sublist), function(db){ 
      
      ps_list[[ds]][[db]] %>% psmelt %>% 
        mutate(PA = if_else(Abundance == 0, 0, 1)) %>% 
        dplyr::select(Sample, all_of(taxRank), PA) %>% 
        filter(PA == 1) %>% 
        as_tibble %>% 
        mutate(database = db,
               dataset = ds)
      
    }) %>% list_rbind
  }) %>% list_rbind
}

PA <- extract_sample_PA(ps_family.ls, 'Family')

# Compute Jaccard by joining tables based on taxa name
# and return a dataframe for every tool comparison
compute_jaccard <- function(df, tool_pair, taxRank) {
  pa1 <- df %>% dplyr::filter(database == tool_pair[1])
  pa2 <- df %>% dplyr::filter(database == tool_pair[2])
  message(paste(tool_pair[1], '&', tool_pair[2]))

  pa_full <- full_join(pa1, pa2, by = c('Sample', taxRank))
  jaccard <- pa_full %>% 
    # exclude species absent from both tools for any sample
    filter(PA.x == 1 | PA.y == 1) %>% 
    # flag common species: 
    mutate(PA.both = ifelse(PA.x == 1 & PA.y == 1, 1, 0)) %>% 
    group_by(Sample) %>% # count common species:
    summarise(PA.both = sum(PA.both, na.rm = TRUE),
              n_tax = n()) %>% 
    # compute jaccard into tibble: 
    mutate(jaccard = PA.both/n_tax, .keep = 'unused')
  
  # for a symmetric heatmap we duplicate the values
  bind_rows(
    jaccard %>% mutate(tool1 = tool_pair[1], tool2 = tool_pair[2]),
    jaccard %>% mutate(tool1 = tool_pair[2], tool2 = tool_pair[1])
  )
}

jaccard_pairwise_df <- PA %>% 
  group_by(dataset) %>%
  group_modify(~ {
    # Create unique tool pairs within the group
    tools <- unique(.x$database)
    tool_pairs <- c(combn(tools, 2, simplify = FALSE))
    # For each instance of that pair, apply the cccvc_compile function
    map_dfr(tool_pairs, function(pair) compute_jaccard(.x, pair, 'Family'))
  }) %>%
  ungroup() %>% # NA on same-tool pairs
  mutate(jaccard = case_when(tool1 == tool2 ~ NA,
                         TRUE ~ jaccard))


