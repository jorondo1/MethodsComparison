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
      dplyr::select(Sample, Grouping_var)
    
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
      pivot_longer(
        cols = c('Richness', 'Shannon', 'Tail', 'Simpson', 'Hill_1', 'Hill_2'),
        names_to = 'Index',
        values_to = 'Index_value') %>% 
      left_join( # Add grouping variable
        samdat, by = 'Sample'
      )
    
  }) %>% list_rbind()
}) %>% list_rbind()

# WILCOX TESTS !
# Set up parallel backend - detect available cores
plan(multisession, workers = availableCores())

alpha_div[['wilcox_tests']] <- alpha_div[['plot_data']] %>%
  group_by(Dataset, Database, Index) %>%
  group_split() %>%
  future_map_dfr(function(df) {
    grp <- df %>% distinct(Dataset, Database, Index)
    
    wt <- wilcox_test(df, as.formula("Index_value ~ Grouping_var"),
                      p.adjust.method = NULL) %>%
      add_significance()
    
    we <- wilcox_effsize(df, as.formula("Index_value ~ Grouping_var"))
    
    bind_cols(grp, wt, we %>% dplyr::select(effsize, magnitude))
  }, .options = furrr_options(seed = TRUE))

# Reset to sequential when done
plan(sequential)

write_rds(alpha_div, 'Out/_Rdata/alpha_div.RDS', compress = 'gz')

# N data naming for datasets (for figures)
dataset_names_n_tibble <- alpha_div[['wilcox_tests']] %>% 
  distinct(Dataset, n) %>% 
  group_by(Dataset) %>% 
  filter(n==max(n)) %>% 
  left_join(tibble(
    Dataset = names(dataset_names),
    Dataset_name = unname(dataset_names)
  ), by = "Dataset") %>% 
  mutate(Dataset_name_n = paste0(Dataset_name, "\n(n=",n,")" ))

# values:
dataset_names_n_lookup <- dataset_names_n_tibble$Dataset_name_n
# add standard names:
names(dataset_names_n_lookup) <- dataset_names_n_tibble$Dataset 
# reorder and flush NAs
dataset_names_n_lookup %<>% .[names(dataset_names)] %>% .[!is.na(.)]
write_rds(dataset_names_n_lookup, 'Out/_Rdata/dataset_names_n_lookup_alphadiv.RDS')


# 2. Beta Diversity ----------------------------------------------------------

# 2.1. Compute bray-curtis PCoA :

pcoa.ls <- list()
future::plan(multicore, workers = detectCores())

pcoa.ls <- map(ps_rare.ls, function(ps.ls) {
  future_map(ps.ls, 
             ~mgx.tools::compute_pcoa(., dist = 'bray', all_coordinates = TRUE), 
             .options = furrr_options(packages = c("phyloseq", "mgx.tools")))
}); plan(sequential)

write_rds(pcoa.ls, 'Out/_Rdata/pcoa_noVST.ls.RDS', compress = 'bz2')

# ompile pairwise distances in a long dataset
pairwise_distances <- imap(pcoa.ls, function(dataset.ls, dataset) {
  imap(dataset.ls, function(database.ls, database) {
    
    compile_dist_pairs(database.ls[['dist.mx']]) %>% 
      mutate(Dataset = dataset,
             Database = database) 
    
  }) %>% list_rbind()
}) %>% list_rbind() 

write_rds(pairwise_distances, 'Out/_Rdata/pairwise_dist.RDS', compress = 'bz2' )


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
      R2 = res$R2[1],   
      p = res$`Pr(>F)`[1]
    ) %>% return()
    
  }) %>% list_rbind()
}) %>% list_rbind()

write_rds(permanova.ds, 'Out/_Rdata/permanova.ds.RDS',compress = 'gz')

