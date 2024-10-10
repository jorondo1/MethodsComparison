# intersect in taxonomy (overlap) and
# pairwise Jaccard index on detected species -> similarity matrix
library(pacman)
p_load(magrittr, tidyverse, phyloseq)

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
                 'NAFLD', 
                 'MPA_db2022', 'MPA_db2023', 
                 Family)


#########################################################
# Pairwise directional taxonomic overlap ###############
# Proportion of taxa found by a tool that are #########
# also found by another tool #########################
#####################################################

taxRank <- 'Genus'

# Compute overlap precision between tool pairs
# Uses taxa_tool_pairs function
# This function applies to a single tool pair subset df
compute_overlap <- function(df, tool_pair, taxRank) {
  pa_full <- df %>% 
    taxa_tool_pairs(tool_pair, taxRank) %>% # custom function
    group_by(Sample) %>% # count common species:
    summarise(PA.xy = sum(Abundance.x*Abundance.y, na.rm = TRUE),
              PA.x = sum(Abundance.x, na.rm = TRUE), # count with tool1
              PA.y = sum(Abundance.y, na.rm = TRUE)) # count with tool2
  
  # return two-way overlap (1 line per comparison)
  bind_rows( 
    pa_full %>% mutate(overlap = PA.xy/PA.x,
                       tool1 = tool_pair[1], 
                       tool2 = tool_pair[2]), 
    pa_full %>% mutate(overlap = PA.xy/PA.y, 
                       tool1 = tool_pair[2],
                       tool2 = tool_pair[1]) 
    # in the resulting df, 'tool1' is the reference taxa set, so the overlap
    # is the proportion of tool1's taxa set also found by tool2
  ) %>% mutate(overlap = case_when(tool1 == tool2 ~ NA, 
                                   TRUE ~ overlap))
}

# Compute overlap for all !
overlap_df <- ps_genus.ls %>% 
  # Melt entire ps list :
  melt_ps_list_glom(taxRank) %>% 
  # Convert abundances to P/A
  mutate(Abundance = if_else(Abundance == 0, 0, 1)) %>% 
  dplyr::select(Sample, all_of(taxRank), Abundance, dataset, database) %>% 
  # iterate overlap calculation over all tool pairs :
  apply_ds_toolpairs(compute_overlap, taxRank) 

# Reorder factors, subset dataset
overlap_df %<>% 
  mutate(dataset = factor(dataset, levels = my_datasets_factorlevels)) 

# plot ! proportion of tool2 (x facet) taxa detected by tool1 (x axis)
overlap_df %>% 
  filter(dataset != 'Moss') %>% 
  ggplot(aes(x = overlap, y = tool1, fill = tool2, colour = tool2)) +
  geom_density_ridges(scale = 0.9, alpha = 0.4, 
                      stat = "binline", boundary = 0, draw_baseline = FALSE
                      ) + xlim(0,1) +
 # facet_grid(~dataset, scales = 'free') +
  ggtitle(paste0('Proportion of taxa identified by other tools in samples (',taxRank,'-level)')) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.15)))+ 
  plot_theme() 

ggsave('Out/overlap_NAFLD_genus.pdf', bg = 'white', 
       width = 2400, height = 1600, units = 'px', dpi = 180)

################
### JACCARD #####
##################

# Jaccard may not be best. Could show proportion of database1
# taxa found by database_[2-n] separately ?

### compute jaccard on a single tool-pair subset df
compute_jaccard <- function(df, tool_pair, taxRank) {
  pa_full <- df %>% 
    taxa_tool_pairs(tool_pair, taxRank) %>% 
    group_by(Sample) %>% # count common species:
    summarise(PA.xy = sum(Abundance.x*Abundance.y, na.rm = TRUE),
              n_tax = n()) %>% 
    # compute jaccard into tibble: 
    mutate(jaccard = PA.xy/n_tax, .keep = 'unused')
  
  # duplicate values to have symmetry (each tool needs to be in position 1 for the plot)
  bind_rows(
    pa_full %>% mutate(tool1 = tool_pair[1], tool2 = tool_pair[2]),
    pa_full %>% mutate(tool1 = tool_pair[2], tool2 = tool_pair[1])
  )
}

jaccard_df <- ps_genus.ls %>% 
  # Melt entire ps list :
  melt_ps_list_glom(taxRank) %>% 
  # Convert abundances to P/A
  mutate(Abundance = if_else(Abundance == 0, 0, 1)) %>% 
  dplyr::select(Sample, all_of(taxRank), Abundance, dataset, database) %>% 
  # iterate overlap calculation over all tool pairs :
  apply_ds_toolpairs(compute_jaccard, taxRank) 

# Reorder factors, subset dataset 
jaccard_df %<>% 
  mutate(dataset = factor(dataset, levels = my_datasets_factorlevels))

# plot ! 
jaccard_df %>% 
  filter(dataset != 'Moss') %>% 
  ggplot(aes(x = jaccard, y = tool1, fill = tool2, colour = tool2)) +
  geom_density_ridges(scale = 0.9, alpha = 0.4, 
                      stat = "binline", boundary = 0, draw_baseline = FALSE
  ) + xlim(0,1)+
  #facet_grid(~dataset, scales = 'free') +
  ggtitle(paste0('Proportion of taxa union set identified by two tools in sample (',taxRank,'-level)')) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.15)))+ 
  plot_theme() 

ggsave('Out/intersect_jaccard_genus.pdf', bg = 'white', 
       width = 2400, height = 1600, units = 'px', dpi = 180)
