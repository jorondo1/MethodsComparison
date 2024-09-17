source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/community_functions.R'))
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/rarefy_even_depth2.R'))

# Bland-Altman analysis, compute the mean and difference between 2 tools
compute_meandiff <- function(div, tool1, tool2) {
  dplyr::select(div, !!sym(tool1), !!sym(tool2)) %>% 
    transmute(mean = (!!sym(tool1) + !!sym(tool2))/2,
              Diff = !!sym(tool1) - !!sym(tool2))
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

# compute multiple diversity indices, output in sublists
div.fun <- function(ps, idx) {
  div_estimate <- list() #initiate list
  for (i in seq_along(idx)) { # compute every diversity index 
    H_q=paste0("H_",i-1) # format H_0, H_1...
    div_estimate[[H_q]] <- estimate_Hill(ps, idx[i])
  }
  return(div_estimate)
}
# compute Hill numbers
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
