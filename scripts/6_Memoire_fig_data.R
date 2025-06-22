library(pacman)
p_load( magrittr, tidyverse, purrr,  phyloseq,
        rstatix, vegan)
ps_rare.ls <- read_rds('Out/ps_rare.ls.rds')
source("scripts/myFunctions.R")
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R'))

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

write_rds(tax_assign_ds, 'Out/_Rdata/tax_assign_ds.RDS')

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

write_rds(tax_assign_sam, 'Out/_Rdata/tax_assign_sam.RDS')


####################
# Alpha diversity ###
######################

# This function will recode them dynamically as A and B, simply 
recode_factor_AB <- function(factor_var) {
  # Ensure variable has exactly 2 levels
  if (nlevels(factor_var) != 2) {
    stop("The input factor must have exactly 2 levels.")
  }
  
  # Recode the factor to "A" and "B" based on level numbers
  factor(as.numeric(factor_var), levels = 1:2, labels = c("A", "B"))
}

alpha_div <- list()
alpha_div[['plot_data']] <-  imap(ps_rare.ls, function(ps_dataset.ls, dataset){
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

# Wilcox test
alpha_div[['wilcox_tests']] <- alpha_div[['plot_data']] %>% 
  group_by(Dataset, Database, Index) %>% 
  wilcox_test(as.formula("Index_value ~ Grouping_var"),
              p.adjust.method = NULL) %>%  # because we want to see what happens when you do only one
  add_significance()

write_rds(alpha_div, 'Out/_Rdata/alpha_div.RDS')


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

write_rds(pcoa.ls, 'Out/_Rdata/pcoa_noVST.ls.RDS')

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
      
    }) %>% list_rbind()
  }) %>% list_rbind()
}) %>% list_rbind() # Add metadata

write_rds(pairwise_distances, 'Out/_Rdata/pairwise_dist.RDS')

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
      ) %>% return()
      
    }) %>% list_rbind()
  }) %>% list_rbind()
}) %>% list_rbind() %>% 
  left_join(CCE_metadata, by = 'Database')

write_rds(permanova.ds, 'Out/_Rdata/permanova.ds.RDS')

