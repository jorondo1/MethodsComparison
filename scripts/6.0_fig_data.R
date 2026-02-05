library(pacman)
p_load( 
  mgx.tools, # devtools::install_github("jorondo1/mgx.tools", force = TRUE)
  magrittr, tidyverse, 
  purrr, parallel, furrr, future,
  phyloseq,
  rstatix, vegan)

ps_rare.ls <- readRDS('Out/_Rdata/ps_rare.ls.rds')
ps_rare.ls$Olive <- NULL
ps_rare.ls %<>% purrr::compact()

source("scripts/myFunctions.R")
source("scripts/0_Config.R")

# Work from the species-level table
taxRanks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")


# 1. Alpha diversity ---------------------------------------------------------

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
               as.factor() %>% recode_factor_AB()) %>% 
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

write_rds(alpha_div, 'Out/_Rdata/alpha_div.RDS', compress = 'gz')


# 2. Beta Diversity ----------------------------------------------------------

# 2.1. Compute bray-curtis PCoA :

pcoa.ls <- list()
future::plan(multicore, workers = detectCores())

pcoa.ls <- map(ps_rare.ls, function(ps.ls) {
  future_map(ps.ls, 
             ~mgx.tools::compute_pcoa(., dist = 'bray', all_coordinates = TRUE), 
             .options = furrr_options(packages = c("phyloseq", "mgx.tools")))
}); plan(sequential)

write_rds(pcoa.ls, 'Out/_Rdata/pcoa_noVST.ls.RDS', compress = 'gz')

# ompile pairwise distances in a long dataset
pairwise_distances <- imap(pcoa.ls, function(dataset.ls, dataset) {
  imap(dataset.ls, function(database.ls, database) {
    
    compile_dist_pairs(database.ls[['dist.mx']]) %>% 
      mutate(Dataset = dataset,
             Database = database) 
    
  }) %>% list_rbind()
}) %>% list_rbind() 

write_rds(pairwise_distances, 'Out/_Rdata/pairwise_dist.RDS', compress = 'gz' )


# 2.2. Procrustes analyses

num_k <- 3

plan(multisession, workers = detectCores())
procrustes_tests <- imap(pcoa.ls, function(dataset.ls, dataset) {
  
  # Create all db combinations
  databases <- names(dataset.ls)
  database_pairs <- expand.grid(databases, databases) %>% 
    filter(Var1 != Var2)
  
  # Parallelized version using future_pmap
  future_pmap(database_pairs, function(Var1, Var2){
    # Extract coordinates
    ds1 <- dataset.ls[[Var1]]$coordinates[,1:num_k]
    ds2 <- dataset.ls[[Var2]]$coordinates[,1:num_k]
    
    # Order / subset samples
    samples_to_keep <- intersect(rownames(ds1), rownames(ds2))
    
    # Procruste test
    test <- protest(ds1[samples_to_keep,], ds2[samples_to_keep,])
    
    # Output 
    tibble(
      ds = dataset,
      db1 = Var1,
      db2 = Var2,
      cor = test$scale,
      pval = test$signif
    )
  }, .options = furrr_options(seed = TRUE)) %>% list_rbind()
  
}) %>% list_rbind(); plan('sequential')

write_rds(procrustes_tests, 'Out/_Rdata/procrustes.RDS',compress = 'gz')

# permanova.ds <- read_rds('Out/permanova.ds.RDS')

# 2.3. perMANOVA
permanova.ds <- imap(pcoa.ls, function(pcoa_db.ls, ds) {
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
}) %>% list_rbind() %>% 
  left_join(CCE_metadata, by = 'Database')

write_rds(permanova.ds, 'Out/_Rdata/permanova.ds.RDS')

