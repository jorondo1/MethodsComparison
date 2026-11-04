# Alpha and beta diversity of ONE dataset -> Out/_Rdata/cache/<DS>/ (alpha_plot_data, alpha_wilcox, pcoa,
# pairwise_dist, procrustes, permanova .RDS). Reads cache/<DS>/ps_rare.RDS (scripts/0_proc_data.R).
# Usually run by scripts/0_run_all.R, which then assembles all datasets into the usual Out/_Rdata/*.RDS files
# (and builds dataset_names_n_lookup_alphadiv.RDS, which needs all datasets).
# By hand: Rscript scripts/1_diversity_metrics.R <dataset>   (or set DS <- '<dataset>' and source() this file)
library(pacman)
p_load( 
  mgx.tools, # devtools::install_github("jorondo1/mgx.tools", force = TRUE)
  magrittr, tidyverse, 
  purrr, parallel, furrr, future,
  phyloseq,
  rstatix, vegan)

if (length(commandArgs(trailingOnly = TRUE)) > 0) DS <- commandArgs(trailingOnly = TRUE)[1]
if (!exists('DS')) stop('No dataset given: Rscript scripts/1_diversity_metrics.R <dataset>')
cache_dir <- file.path('Out/_Rdata/cache', DS)
message('### ', DS, ': alpha and beta diversity -> ', cache_dir)

# A one-dataset list, so that the imap(dataset) -> imap(database) code below runs as it did on all datasets
ps_rare.ls <- list(readRDS(file.path(cache_dir, 'ps_rare.RDS'))) %>% set_names(DS)
ps_rare.ls[[DS]] %<>% purrr::compact()

source("scripts/0_config.R")

# Datasets without a grouping variable get alpha diversity, PCoA, pairwise distances and Procrustes, but no Wilcoxon
# tests or PERMANOVA. Samples with no group (e.g. the Maple samples without Provenance) keep their diversity metrics
# and are only left out of the tests.
has_groups <- function(ds) ds %in% names(grouping_variable) && !is.na(grouping_variable[[ds]])

# Recode a two-level factor as A and B (was in scripts/archive/myFunctions.R; everything else comes from mgx.tools)
recode_factor_AB <- function(factor_var) {
  if (nlevels(factor_var) != 2) stop("The input factor must have exactly 2 levels.")
  factor(as.numeric(factor_var), levels = 1:2, labels = c("A", "B"))
}

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
      { if (has_groups(dataset)) {
          mutate(., Grouping_var = !!sym(grouping_variable[[dataset]]) %>% as.factor() %>% recode_factor_AB())
        } else mutate(., Grouping_var = factor(NA, levels = c('A', 'B'))) } %>% 
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
# Set up parallel backend: availableCores() = the Slurm allocation on the HPC (detectCores() would count the whole node)
plan(multisession, workers = availableCores())

alpha_div[['wilcox_tests']] <- alpha_div[['plot_data']] %>%
  filter(map_lgl(Dataset, has_groups), !is.na(Grouping_var)) %>%
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

write_rds(alpha_div[['plot_data']], file.path(cache_dir, 'alpha_plot_data.RDS'), compress = 'gz')
write_rds(alpha_div[['wilcox_tests']], file.path(cache_dir, 'alpha_wilcox.RDS'), compress = 'gz')

# (Dataset labels with n for the figures: built from all datasets in scripts/0_run_all.R)


# 2. Beta Diversity ----------------------------------------------------------

# 2.1. Compute bray-curtis PCoA :

pcoa.ls <- list()
future::plan(multicore, workers = availableCores())

pcoa.ls <- map(ps_rare.ls, function(ps.ls) {
  future_map(ps.ls, 
             ~mgx.tools::compute_pcoa(., dist = 'bray', all_coordinates = TRUE), 
             .options = furrr_options(packages = c("phyloseq", "mgx.tools")))
}); plan(sequential)

write_rds(pcoa.ls[[DS]], file.path(cache_dir, 'pcoa.RDS'), compress = 'gz')

# ompile pairwise distances in a long dataset
pairwise_distances <- imap(pcoa.ls, function(dataset.ls, dataset) {
  imap(dataset.ls, function(database.ls, database) {
    
    compile_dist_pairs(database.ls[['dist.mx']]) %>% 
      mutate(Dataset = dataset,
             Database = database) 
    
  }) %>% list_rbind()
}) %>% list_rbind() 

write_rds(pairwise_distances, file.path(cache_dir, 'pairwise_dist.RDS'), compress = 'gz')


# 2.2. Procrustes analyses

num_k <- 3

plan(multisession, workers = availableCores())
set.seed(1234) # furrr derives each pair's permutation stream from it: reproducible p-values
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

write_rds(procrustes_tests, file.path(cache_dir, 'procrustes.RDS'), compress = 'gz')

# 2.3. perMANOVA
permanova.ds <- imap(pcoa.ls[map_lgl(names(pcoa.ls), has_groups)], function(pcoa_db.ls, ds) {
  imap(pcoa_db.ls, function(pcoa, db) { # iterate over distances
    
    samData <- pcoa$metadata %>% 
      mutate(Grouping_var = factor(!!sym(grouping_variable[[ds]]), 
                                   labels = c("Group 1", "Group 2")))
    # samples without a group are left out of the test (metadata rows are in the distance matrix's order)
    grouped <- !is.na(samData$Grouping_var)
    dist.mx <- as.dist(as.matrix(pcoa$dist.mx)[grouped, grouped])
    samData <- samData[grouped, , drop = FALSE]
    
    # permanova (seeded: the permutations, hence p, are reproducible)
    set.seed(1234)
    res <- adonis2(formula = dist.mx ~ Grouping_var, 
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

write_rds(permanova.ds, file.path(cache_dir, 'permanova.RDS'), compress = 'gz')

