library(pacman)
p_load(magrittr, tidyverse, phyloseq)

source('scripts/myFunctions.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/tax_glom2.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/psflashmelt.R')

#########################################################
### Scaled difference in mean group relative abundance ###
###########################################################

ps_filt.ls <- readRDS('Out/_Rdata/ps_filt.ls.RDS')

which_taxrank <- c('Genus', 'Family')

mean_taxa_ratios_groups <- function(ps, which_ds, taxRank) {
  relAb_table <- ps %>% 
    tax_glom2(taxrank = taxRank) %>% 
    psflashmelt() %>% 
    #    filter(Abundance>0) %>%  # Don't !!
    group_by(Sample) %>% 
    dplyr::mutate(relAb = Abundance / sum(Abundance),
                  Taxon = OTU) %>% 
    ungroup() %>% # otherwise the next mutate reinstates the Sample variable in the grouping_variable vector
    # Recode grouping as AB
    mutate(Grouping_var = !!sym(grouping_variable[which_ds]) %>% 
             as.factor() %>% recode_factor_AB()
           ) %>% 
    dplyr::select(Sample, Taxon, relAb, Grouping_var)
  
  # Compute relative abundance score (mean relative abundance of each taxon across samples)
  relAb_scores <- relAb_table %>% 
    group_by(Taxon) %>% 
    dplyr::summarise(meanRelAb = mean(relAb))
  
  # Compute scaled differences
  scaled_diffs <- relAb_table %>% 
    group_by(Grouping_var, Taxon) %>% 
    dplyr::summarise(meanRelAb = mean(relAb),
                     sdRelAb = sd(relAb),
                     prev = sum(relAb != 0),
                     .groups = 'drop') %>% 
    group_by(Taxon) %>% 
    # Ratios get tricky with low abundance taxa because quasi absence from a group inflates the ratios drastically (even log)
    dplyr::summarise(scaled_diff =
                       (max(meanRelAb[Grouping_var == "A"],
                            meanRelAb[Grouping_var == "B"]) -
                          min(meanRelAb[Grouping_var == "A"],
                              meanRelAb[Grouping_var == "B"])) /
                       max(meanRelAb[Grouping_var == "A"],
                           meanRelAb[Grouping_var == "B"])
    )
  
  # One output joined
  left_join(scaled_diffs, relAb_scores, by = 'Taxon') %>% 
    return()
}

out_list <- list()

for (taxRank in which_taxrank) {
  out_list[[taxRank]] <- 
    imap(ps_filt.ls, function(ps_list, dataset) {
      imap(ps_list, function(ps, database) {
        mean_taxa_ratios_groups(ps, dataset, taxRank) %>% 
          mutate(
            database = database
          )
      }) %>% list_rbind()
    })
}

write_rds(out_list, 'Out/_Rdata/taxa_relAb_metrics.RDS')
