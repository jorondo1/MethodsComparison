# intersect in taxonomy
# pairwise Jaccard index on detected species -> similarity matrix
# First test at genus level
library(pacman)
p_load(magrittr, tidyverse, phyloseq, rlang)

ps_genus.ls <- read_rds("Out/ps_rare_genus.ls.rds") 
ps_family.ls <- read_rds("Out/ps_rare_family.ls.rds") 

##############
# Intersect ###
##############

# return taxa from exclusion space
extract_taxnames <- function(ps, ds, db1, db2, taxRank) {
  
   extract_tax <- function(db) {
     ps[[ds]][[db]] %>% psmelt %>% 
       pull({{ taxRank }}) %>% unique # extract unique ranks vector
   }
  tax1 <- extract_tax(db1)
  tax2 <- extract_tax(db2)
  intersect_set <- intersect(tax1,tax2)
  message(paste(length(intersect_set), 'found in common'))
  union_set <- union(tax1, tax2) 
  message(paste(length(union_set), 'found in total'))
  setdiff(union_set, intersect_set) %>% sort
}

extract_taxnames(ps_genus.ls,
                 'Saliva', 
                 'SM_genbank_202203', 'SM_gtdb_rs214_full', 
                 Genus)

taxRank <- 'Family'

# Create a Presence/Absence matrix for all datasets & tools
extract_sample_PA <- function(ps_list, taxRank) {
  map(names(ps_list), function(ds){ # first list level: apply to all dataset
    ds_sublist <- ps_list[[ds]] # extract sublist (ps object)
    map(names(ds_sublist), function(db){ 
      # summarize by higher level taxonomy (not always useful, having a pre-summarized phyloseq is better)
      ps_list[[ds]][[db]] %>% psmelt %>% 
        group_by(Sample, !!sym(taxRank)) %>% 
        summarise(Abundance = sum(Abundance), 
                  .groups = 'drop') %>% 
        # Convert to P/A
        mutate(PA = if_else(Abundance == 0, 0, 1)) %>% 
        dplyr::select(Sample, all_of(taxRank), PA) %>% 
        filter(PA == 1) %>% # reduce matrix size; absences can be removed
        as_tibble %>% 
        mutate(database = db,
               dataset = ds)
    }) %>% list_rbind
  }) %>% list_rbind
}

PA <- extract_sample_PA(ps_family.ls, taxRank)

#########################################################
# Pairwise directional taxonomic overlap ###############
# Proportion of taxa found by a tool that are #########
# also found by another tool #########################
#####################################################

overlap_precision <- function(df, tool_pair, taxRank) {
  pa_full <- df %>% 
    taxa_tool_pairs(tool_pair, taxRank) %>% # custom function
    group_by(Sample) %>% # count common species:
    summarise(PA.both = sum(PA.both, na.rm = TRUE),
              PA.x = sum(PA.x, na.rm = TRUE), # count with tool1
              PA.y = sum(PA.y, na.rm = TRUE)) # count with tool2
  
  # return two-way overlap (1 line per comparison)
  bind_rows( 
     pa_full %>% mutate(overlap = PA.both/PA.x,
                        tool1 = tool_pair[1], 
                        tool2 = tool_pair[2]), 
     pa_full %>% mutate(overlap = PA.both/PA.y, 
                          tool1 = tool_pair[2],
                          tool2 = tool_pair[1])
     # in the resulting df, 'tool1' is the reference taxa set, so the overlap
     # is the proportion of tool1's taxa set also found by tool2
  )
}

# Apply to each dataset and iterate over tool pairs
overlap_df <- PA %>% 
  group_by(dataset) %>%
  group_modify(~ {
    # Create unique tool pairs within the group
    tools <- unique(.x$database)
    tool_pairs <- c(combn(tools, 2, simplify = FALSE))
    # For each instance of that pair, apply the cccvc_compile function
    map_dfr(tool_pairs, function(pair) overlap_precision(.x, pair, taxRank))
  }) %>%
  ungroup() %>% # NA on same-tool pairs
  mutate(overlap = case_when(tool1 == tool2 ~ NA,
                         TRUE ~ overlap))

# Reorder factors, subset dataset
overlap_df %<>% 
  mutate(dataset = factor(dataset, levels = c('Saliva', 'Feces', 'Moss'))) %>% 
  filter(dataset != 'Moss')

# plot ! proportion of tool2 (x facet) taxa detected by tool1 (x axis)
overlap_df %>% 
  ggplot(aes(x = tool2, y = overlap, fill = tool2)) +
  geom_violin() +
  facet_grid(dataset ~ tool1, scales = 'free_x') +
  ylim(0,1) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle(paste('Proportion of taxa identified by other tools in samples (',taxRank,'-level)'))

ggsave('Out/intersect_overlap.pdf', bg = 'white', 
       width = 2400, height = 1600, units = 'px', dpi = 180)

# Jaccard may not be best. Could show proportion of database1
# taxa found by database_[2-n] separately ?

### JACCARD
compute_jaccard <- function(df, tool_pair, taxRank) {
  pa_full <- df %>% 
    taxa_tool_pairs(tool_pair, taxRank) %>% 
    group_by(Sample) %>% # count common species:
    summarise(PA.both = sum(PA.both, na.rm = TRUE),
              n_tax = n()) %>% 
    # compute jaccard into tibble: 
    mutate(jaccard = PA.both/n_tax, .keep = 'unused')
  
  # duplicate values to have symmetry (each tool needs to be in position 1 for the plot)
  bind_rows(
    pa_full %>% mutate(tool1 = tool_pair[1], tool2 = tool_pair[2]),
    pa_full %>% mutate(tool1 = tool_pair[2], tool2 = tool_pair[1])
  )
}

# Compute pairwise jaccard across all 
jaccard_pairwise_df <- PA %>% 
  group_by(dataset) %>%
  group_modify(~ {
    # Create unique tool pairs within the group
    tools <- unique(.x$database)
    tool_pairs <- c(combn(tools, 2, simplify = FALSE))
    # For each instance of that pair, apply the cccvc_compile function
    map_dfr(tool_pairs, function(pair) compute_jaccard(.x, pair, taxRank))
  }) %>%
  ungroup() %>% # NA on same-tool pairs
  mutate(jaccard = case_when(tool1 == tool2 ~ NA,
                             TRUE ~ jaccard))

# Reorder factors, subset dataset
jaccard_pairwise_df %<>% 
  mutate(dataset = factor(dataset, levels = c('Saliva', 'Feces', 'Moss'))) %>% 
  filter(dataset != 'Moss')

# plot ! 
jaccard_pairwise_df %>% 
  ggplot(aes(x = tool1, y = jaccard, fill = tool1)) +
  geom_violin() +
  facet_grid(dataset ~ tool2, scales = 'free_x') +
  ylim(0,1) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

ggsave('Out/intersect_jaccard.pdf', bg = 'white', 
       width = 2400, height = 1600, units = 'px', dpi = 180)
