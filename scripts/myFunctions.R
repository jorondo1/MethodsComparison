source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/community_functions.R'))
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/rarefy_even_depth2.R'))

filter_low_prevalence <- function(ps, minPrev = 0.05, minAbund = 0.001) {
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
    summarise(across(where(is.numeric), sum)) %>% 
    column_to_rownames('Species') 
  
  # Extract taxonomy
  tax <- abunTable %>% 
    dplyr::select(where(is.character)) %>% 
    unique %>% # because of renaming above, some species will be duplicate
    mutate(Species2 = Species) %>% 
    column_to_rownames('Species2') %>% as.matrix
  
  # Build phyloseq
  phyloseq(otu_table(abund, taxa_are_rows = TRUE),
           sample_data(sampleData),
           tax_table(tax)
  ) %>% 
    (if (filtering) filter_low_prevalence else identity)
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

# Hill numbers
estimate_Hill <- function(ps, q) {
  x <- ps@otu_table %>% as("matrix")
  if (taxa_are_rows(ps)) { 
    x <- t(x) 
  }
  total <- rowSums(x)
  x <- sweep(x, 1, total, "/")
  
  if (q == 0) {  # Species richness
    div <- rowSums(x > 0)
  } else if (q == 1) { # Shannon diversity (exponential of Shannon entropy)
    div <- exp(-rowSums(x * log(x, base = exp(1)), na.rm = TRUE))
  } else {  # Hill number formula for q ≠ 0 and q ≠ 1
    div <- rowSums(x^q)^(1 / (1 - q))
  }
  return(div)
}

# Esitmate diversity (Shannon, Simpson, Tail)
estimate_diversity <- function(ps, index = 'Shannon') {
  x <- ps@otu_table %>% as("matrix")
  if (taxa_are_rows(ps)) { 
    x <- t(x) 
  }
  total <- apply(x, 1, sum)
  x <- sweep(x, 1, total, "/")
  
  if(index == 'Tail') {
    tail_stat <- function(row) {
      values <- sort(row, decreasing = TRUE)
      sqrt(sum(values * ((seq_along(values)-1)^2)))
    }
    div <- apply(x, 1, tail_stat)
  }
  if(index == 'Shannon') {
    x <- -x * log(x, exp(1))
    div <- apply(x, 1, sum, na.rm = TRUE)
  }
  if(index == 'Simpson') {
    div <- 1 - apply((x * x), 1, sum, na.rm = TRUE) 
  }
  if(index == 'Richness') {
    div <- apply(x, 1, function(x) sum(x != 0))
  }
  return(div)
}

# compute multiple diversity indices, output in sublists
div.fun <- function(ps, idx) {
  div_estimate <- list() #initiate list
  for (i in seq_along(idx)) { # compute Hill numbers
    H_q=paste0("H_",i-1) # format H_0, H_1...
    div_estimate[[H_q]] <- estimate_Hill(ps, idx[i])
  }
  div_estimate[["Tail"]] <- estimate_diversity(ps, index = "Tail")
  return(div_estimate)
}

################
### Other functions 
################

# Filtering dataset for a tool pair and 
# generate a list of taxa found by either tool
taxa_tool_pairs <- function(df, tool_pair, taxRank) {
  pa1 <- df %>% dplyr::filter(database == tool_pair[1])
  pa2 <- df %>% dplyr::filter(database == tool_pair[2])
  
  message(paste(tool_pair[1], '&', tool_pair[2]))
  
  # return the full set of taxa identified by either tool
  full_join(pa1, pa2, by = c('Sample', taxRank)) %>% 
    # exclude species absent from both tools for any sample
    filter(PA.x == 1 | PA.y == 1) %>% 
    # flag common species: 
    mutate(PA.both = PA.x*PA.y)
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

# Bland-Altman analysis, compute the mean and difference between 2 tools
compute_meandiff <- function(div, tool1, tool2) {
  dplyr::select(div, !!sym(tool1), !!sym(tool2)) %>% 
    transmute(mean = (!!sym(tool1) + !!sym(tool2))/2,
              Diff = !!sym(tool1) - !!sym(tool2))
}

