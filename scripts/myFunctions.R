plot_theme <- function() {
  list(
    theme_minimal(),
    scale_fill_manual(values = tool_colours),
    scale_colour_manual(values = tool_colours)
  )
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

# Recode a two-level factor as A and B
recode_factor_AB <- function(factor_var) {
  # Ensure variable has exactly 2 levels
  if (nlevels(factor_var) != 2) {
    stop("The input factor must have exactly 2 levels.")
  }
  
  # Recode the factor to "A" and "B" based on level numbers
  factor(as.numeric(factor_var), levels = 1:2, labels = c("A", "B"))
}


# Melt ps object list, agglomerate abundances to desired taxRank
melt_ps_list_glom <- function(ps_list, taxRank) {
  map(names(ps_list), function(ds){ # first list level: apply to all dataset
    ds_sublist <- ps_list[[ds]] # extract sublist (ps object)
    map(names(ds_sublist), function(db){ 
      # summarize by higher level taxonomy (not always useful, having a pre-summarized phyloseq is better)
      ps_list[[ds]][[db]] %>% psflashmelt %>% 
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

compute_3_lvl <- function(ps.ls, func, ...){
  require('furrr') 
   # not sure if multicore is better 
  
  imap(ps.ls, function(taxRank.ls, taxRank) {
    cat("Processing", taxRank, "...\n")
   # if (taxRank != "Species") return(NULL)      # ! DEV !
    
    imap(taxRank.ls, function(ds.ls, ds) {
      samVar <- grouping_variable[[ds]]               # Group variable to test 
     # if (ds != "NAFLD") return(NULL)          # ! DEV !
      cat("Processing dataset:", ds, "...\n")
      
      # ! Warning : doesn't work from RStudio! 
      future_imap(ds.ls, function(db.ps, db) { 
        message("Using database:", db, "...\n")
        
        # Collect all available arguments
        all_args <- list(ps = db.ps, 
                         samVar = samVar, 
                         taxRank = taxRank, 
                         ds = ds, 
                         db = db,
                         out_path = out_path)
        
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