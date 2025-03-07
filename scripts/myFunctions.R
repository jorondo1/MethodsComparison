source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/rarefy_even_depth2.R'))
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/phyloseq_to_edgeR.R'))

my_datasets_factorlevels <- c('P19_Saliva', 'P19_Gut', 'RA_Gut', 'AD_Skin', 'Moss', 'NAFLD')
tool_colours <- c(
  'MPA_db2019' = 'seagreen3',
  'MPA_db2022' = 'darkolivegreen2',
  'MPA_db2023' = 'darkolivegreen4',
  'MOTUS' = 'goldenrod',
  'KB10' = 'indianred1',
  'KB45' = 'indianred2',
  'KB51' = 'indianred3',
  'KB90' = 'orangered4',
  'SM_genbank-2022.03' = 'purple3',
  'SM_gtdb-rs214-full' = 'navyblue',
  'SM_gtdb-rs214-rep'= 'royalblue',
  'SM_gtdb-rs214-rep_MAGs'= 'skyblue3'
)

group_vars <- c(
  'NAFLD' = 'Group',
  'AD_Skin' = 'Group',
  'Moss' = 'Compartment'
)

CCE_names <- c(
  'MOTUS' = 'mOTUs',
  'MPA_db2019' = 'MetaPhlAn3 (2019)',
  'MPA_db2022' = 'Metaphlan4 (2022)',
  'MPA_db2023' = 'Metaphlan4 (2023)',
  'KB10' = 'Kraken (10% conf.)\n+ Bracken',
  'KB45' = 'Kraken (45% conf.)\n+ Bracken',
  'KB51' = 'Kraken (51% conf.)\n+ Bracken',
  'KB90' = 'Kraken (90% conf.)\n+ Bracken',
  'SM_genbank-2022.03' = 'Sourmash (Genbank)',
  'SM_gtdb-rs214-full' = 'Sourmash (GTDB)',
  'SM_gtdb-rs214-rep'= 'Sourmash (GTDB rep.)',
  'SM_gtdb-rs214-rep_MAGs'= 'Sourmash (GTDB rep.\n+ Novel MAGs)'
)

Hill_numbers <- c(
  'H_0' = 'Richness',
  'H_1' = 'Shannon (Hill order 1)',
  'H_2' = 'Simpson (Hill order 2)'
)

CCE_metadata <- tibble(
  database = c("MOTUS","MPA_db2019","MPA_db2022","MPA_db2023",
               "KB10","KB45","KB51","KB90",
               "SM_genbank-2022.03", "SM_gtdb-rs214-full","SM_gtdb-rs214-rep","SM_gtdb-rs214-rep_MAGs"),
  plot_colour = c("goldenrod","green4","palegreen3","darkolivegreen",
                  "indianred1","indianred2","indianred3", "orangered4", 
                  "purple3", "navyblue", "royalblue", "skyblue3"),
  tool = c('MOTUS', 'MPA','MPA','MPA',
           'KB','KB','KB','KB','SM','SM','SM','SM'),
  CCE_approach = c('DNA-to-Marker','DNA-to-Marker','DNA-to-Marker','DNA-to-Marker',
               'DNA-to-DNA','DNA-to-DNA','DNA-to-DNA','DNA-to-DNA','DNA-to-DNA','DNA-to-DNA','DNA-to-DNA','DNA-to-DNA'),
  taxonomy = c('NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','GTDB','GTDB','GTDB'),
  
)


tool_vars <- tibble(
  "Aldex2" = "wi.eBH",
  "ANCOMBC2" = "adj_p",
  "radEmu" = "pval",
  "corncob" = "p_fdr",
  "DESeq2" = "padj",
  "edgeR" = "FDR", 
  "MaAsLin2" = "qval"
)

DAA_metadata <- tibble(
  DAA_tool = names(tool_vars),
  Compositional = c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE),
  Distribution = c('DM', 'LL', 'LL', 'BB', 'NB', 'NB', 'N'),
  Transformation = c('CLR', 'LT', 'LT', 'NONE', 'NONE', 'NONE', 'AST'),
  Taxon_bias = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  Target_ = c('RA', 'AA','AA','RA','RA','RA','RA'), # Absolute or relative abundance
  plot_shape = c(15, 19, 17, 18, 4, 3, 6)
)


plot_theme <- function() {
  list(
    theme_minimal(),
    scale_fill_manual(values = tool_colours),
    scale_colour_manual(values = tool_colours)
  )
}

filter_low_prevalence <- function(ps, minPrev = 0.05, minAbund = 0.0001) {
  # Taxa prevalence
  prev <- apply(otu_table(ps), 1, function(x) sum(x > 0)) / nsamples(ps)
  
  # Convert to relative abundance
  rel_abund <- apply(otu_table(ps), 2, function(x) x / sum(x))
  
  # Keep taxa with prevalence above threshold
  keepTaxa <- names(prev[prev >= minPrev & apply(rel_abund, 1, max) >= minAbund])
  
  # Subset phyloseq object
  prune_taxa(keepTaxa, ps) %>% return
}

### Build phyloseq object from MPA output
assemble_phyloseq <- function(abunTable, sampleData, filtering = FALSE) {
  
  abunTable %<>% dplyr::filter(Kingdom == "Bacteria") %>% 
    mutate(across(where(is.character), \(x) {
      str_replace_all(x,'_', ' ') %>%
        str_replace('Candidatus ', '') %>% 
        str_remove(" [A-Z]$")  #https://gtdb.ecogenomic.org/faq#why-do-some-family-and-higher-rank-names-end-with-an-alphabetic-suffix
    }))
    
  # Extract abundance table with Species as identifier
  abund <- abunTable %>% 
    dplyr::select(where(is.double), Species) %>% 
    group_by(Species) %>% 
    dplyr::summarise(across(where(is.numeric), sum)) %>% 
    column_to_rownames('Species') %>% 
    select(where(~ sum(.) >= 100))
  
  # Extract taxonomy
  tax <- abunTable %>% 
    dplyr::select(where(is.character)) %>% 
    unique %>% # because of renaming above, some species will be duplicate
    mutate(Species2 = Species) %>% 
    column_to_rownames('Species2') %>% as.matrix
  
  # Some datasets may end up with very low read counts and lose samples.
  # We subset the sample dataset, but we add a check if all samples are lost:
  keep_samples <- sum(rownames(sampleData) %in% colnames(abund))
  
  if (keep_samples==0) {
    return(NULL)
    } else {
    sampleData_subset <- sampleData[keep_samples,]
    
    # Build phyloseq
    ps <- phyloseq(otu_table(abund, taxa_are_rows = TRUE),
             sample_data(sampleData_subset),
             tax_table(tax)
    ) %>% 
      (if (filtering) filter_low_prevalence else identity)
    
    prune_samples(sample_sums(ps) > 0, ps) %>%  #remove any empty samples 
      prune_taxa(taxa_sums(.) > 0,.) # remove taxa absent from all (may happen if you end up using not all the samples you parse, e.g. metadata missing so sample dropped in the process)
  }
}

# Compute sparseness (proportion of 0 in abundance matrix)
compile_sparseness <- function(ps_list) {
  map(names(ps_list), function(ds){ # first list level: apply to all dataset
    ds_sublist <- ps_list[[ds]] # extract sublist (ps object)
    map(names(ds_sublist), function(db){ # second list : apply to all ps
      ps <- ds_sublist[[db]]
      otu <- otu_table(ps)
      tibble(
        sparseness = sum(otu == 0)/length(otu), # sparseness
        dataset = ds,
        database = db
      )
    }) %>% list_rbind
  }) %>% list_rbind
}


################
### DIVERSITY ###
################

# Add sample data to a dataframe from phyloseq objeoct
filter_and_add_samData <- function(df, ds, Rank, ps_list) {
  df %<>% filter(dataset == !!ds & Rank == !!Rank)
  
  # Extract sam_data from any phyloseq object of that dataset
  samData <- ps_list[[Rank]][[ds]][[1]]@sam_data %>% 
    as('data.frame') %>% 
    rownames_to_column('Sample')
  
  inner_join(df, samData, by = 'Sample') 
}

################
### General functions 
################

# Melt ps object list, agglomerate abundances to desired taxRank
melt_ps_list_glom <- function(ps_list, taxRank) {
  map(names(ps_list), function(ds){ # first list level: apply to all dataset
    ds_sublist <- ps_list[[ds]] # extract sublist (ps object)
    map(names(ds_sublist), function(db){ 
      # summarize by higher level taxonomy (not always useful, having a pre-summarized phyloseq is better)
      ps_list[[ds]][[db]] %>% psmelt %>% 
        group_by(Sample, !!sym(taxRank)) %>% 
        dplyr::summarise(Abundance = sum(Abundance), 
                  .groups = 'drop') %>% 
        dplyr::filter(Abundance !=0) %>% 
        as_tibble %>% 
        mutate(database = db,
               dataset = ds)
    }) %>% list_rbind
  }) %>% list_rbind
}

# Applies a function to each possible combination of tools (database) 
# found in the dataset; passes the pair names to the function
apply_ds_toolpairs <- function(df, func, taxRank) {
  df %>% group_by(dataset) %>%
    group_modify(~ {
      data_subset <- .x
      # Create unique tool pairs within the group
      tools <- unique(data_subset$database)
      tool_pairs <- combn(tools, 2, simplify = FALSE)
      
      # For each instance of that pair, apply the cccvc_compile function
      map_dfr(tool_pairs, function(pair) func(data_subset, pair, taxRank))
    })
}

# Filtering dataset for a tool pair and 
# generate a list of taxa found by either tool
taxa_tool_pairs <- function(df, tool_pair, taxRank) {
  set.x <- df %>% dplyr::filter(database == tool_pair[1])
  set.y <- df %>% dplyr::filter(database == tool_pair[2])
  message(paste(tool_pair[1], '&', tool_pair[2]))
  
  # return the full set of taxa identified by either tool
  full_join(set.x, set.y, by = c('Sample', taxRank)) %>% 
    # exclude species absent from both tools for any sample
    filter(Abundance.x != 0 | Abundance.y != 0) 
}

# Function to apply a function iteratively to every last object in the list,
# and pass it the level names as arguments (db, ds, dist) : 
iterate_distances <- function(pcoa.ls, compile_func) {
  map(names(pcoa.ls), function(ds) { #iterate over dataset names
    ds_sublist <- pcoa.ls[[ds]] 
    map(names(ds_sublist), function(db) { # iterate over databases
      ps <- ds_sublist[[db]]
      map(names(ps), function(dist) { # iterate over distances
        compile_func(ps, ds, db, dist)
      }) %>% list_rbind
    }) %>% list_rbind
  }) %>% list_rbind 
}

# Extract name of lowest available taxRank in a ps object
extract_lowest_rank <- function(ps) {
  require('phyloseq')
  require('dplyr')
  tax_table(ps) %>%
    as.data.frame %>%
    dplyr::summarise(across(everything(), ~ !all(is.na(.)))) %>%
    tidyr::pivot_longer(everything()) %>%
    filter(value) %>%
    slice_tail(n = 1) %>%
    pull(name)
}

# Using imap to manage list names in 3-level lists-of-lists
# 1. taxonomic rank
# 2. dataset (ds)
# 3. tool/database (db)
# Generates a list with the same hierarchy with what func() returns as 
# the lowest-level objects
plan(multicore, workers = 9)
compute_3_lvl <- function(ps.ls, func, ...){
  require('furrr') 
   # not sure if multicore is better 
  
  imap(ps.ls, function(taxRank.ls, taxRank) {
    cat("Processing", taxRank, "...\n")
   # if (taxRank != "Species") return(NULL)      # ! DEV !
    
    imap(taxRank.ls, function(ds.ls, ds) {
      samVar <- group_vars[[ds]]               # Group variable to test 
     # if (ds != "NAFLD") return(NULL)          # ! DEV !
      cat("Processing dataset:", ds, "...\n")
      
      future_imap(ds.ls, function(db.ps, db) { # ! Warning : doesn't work from RStudio! 
        message("Using database:", db, "...\n")
        
        # Collect all available arguments
        all_args <- list(ps = db.ps, samVar = samVar, 
                         taxRank = taxRank, ds = ds, db = db, out_path = out_path)
        # keep required arguments only
        func_args <- all_args[names(all_args) %in% names(formals(func))]
        # Call computing function: 
        do.call(func, func_args)
        
      }, .options = furrr_options(seed = T))
    })
  })
}

compile_3_lvl <- function(results.ls, func, ...) {
  imap(results.ls, function(taxRank.ls, taxRank) {
    imap(taxRank.ls, function(ds.ls, ds) {
      imap(ds.ls, function(db.ps, db) { 
        func(db.ps, taxRank, db, ds) 
      }) %>% list_rbind
    }) %>% list_rbind
  }) %>% list_rbind
}

################
### METRICS ###
################

# Function to apply cccvc to a pair of tools, with error handling
cccvc_compile <- function(df, tool_pair) {
  # Filter the data to keep only the two tools in the pair
  df_pair <- df %>% filter(database %in% tool_pair)
  
  # Try to compute cccvc, and return NA in case of an error
  tryCatch({
    # Run cccvc and return results in a tibble
    ccc_out <- cccvc(df_pair, ry = 'value', rind = 'Sample', rmet = 'database')
    
    tibble(
      tool1 = c(tool_pair[1], tool_pair[2]),
      tool2 = c(tool_pair[2], tool_pair[1]),
      CCC = c(ccc_out$ccc[1], ccc_out$ccc[1]),       # CCC 
      LL_CI_95 = c(ccc_out$ccc[2], ccc_out$ccc[2]),  # Lower limit CI
      UL_CI_95 = c(ccc_out$ccc[3], ccc_out$ccc[3]),  # Upper limit CI
      SE_CCC = c(ccc_out$ccc[4], ccc_out$ccc[4])     # Standard error
      
    )
  }, error = function(e) {
    tibble( # If error, don't store any values except tool names
      tool1 = tool_pair[1],
      tool2 = tool_pair[2]
    )
  })
}

# Pairwise distances
compute_pairwise_dist <- function(sample_df) {
  # Select only numeric columns (tool abundances) to calculate distances
  tool_est.mx <- sample_df %>% select(where(is.numeric)) %>% as.matrix 
  
  # Compute pairwise distances for each combination of tools
  dist_vectors <- list()
  tool_pairs <- combn(ncol(tool_est.mx), 2) # Generate all pairs of columns (tools)
  for (i in 1:ncol(tool_pairs)) {
    tool1 <- tool_pairs[1, i]
    tool2 <- tool_pairs[2, i]
    
    # Extract tool pair as matrix
    tool_pair_matrix <- tool_est.mx[, c(tool1, tool2), drop = FALSE] %>% 
      .[rowSums(.) != 0, ] %>% t
    
    # Compute robust Aitchison distance tools
    dist_value <- vegdist(tool_pair_matrix, method = "robust.aitchison")
    
    # Store distance as vector
    dist_vectors[[i]] <- setNames(as.vector(dist_value), 
                                  paste(colnames(tool_est.mx)[tool1], 
                                        colnames(tool_est.mx)[tool2], 
                                        sep = "_vs_"))
  }
  # Flatten vector list into single vector
  return(unlist(dist_vectors))
} 

# Bland-Altman analysis, compute the mean and difference between 2 tools
compute_meandiff <- function(div, tool1, tool2) {
  dplyr::select(div, !!sym(tool1), !!sym(tool2)) %>% 
    transmute(mean = (!!sym(tool1) + !!sym(tool2))/2,
              Diff = !!sym(tool1) - !!sym(tool2))
}

######################################
### Differential Abundance Testing ####
########################################