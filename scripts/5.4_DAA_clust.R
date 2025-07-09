library(pacman)
p_load(tidyverse, magrittr, ape, RColorBrewer, ComplexHeatmap, wCorr)
source("scripts/myFunctions.R")

# Import data
DAA <- Sys.glob('Out/DAA/*.tsv') %>% 
  map(read_tsv, show_col_types = FALSE) %>% 
  list_rbind()

# Define sets to work with 
which_tools <- c('DESeq2', 'edgeR', 'Aldex2','corncob',
                 'MaAsLin2',  'ZicoSeq', 'ANCOMBC2'#, 'radEMU'
)
which_databases <- c('KB45', 'KB90', 'SM_RefSeq_20250528', 
                     'KB45_GTDB', 'KB90_GTDB', 'SM_gtdb-rs220-rep', 
                     'MPA_db2023', 'MOTUS')


# Function to compute symmetric signed ranks of coefficients (effect sizes)
# across DAA for taxon found in at least <at_least_DAA> methods. 

rank_correlation_by_ds_db <- function(DAA_tibble, at_least_DAA) {
  
  # Keep taxa detected by most methods
  which_taxa <- DAA_tibble %>% 
    group_by(Taxon) %>% 
    summarise(n = n()) %>% 
    filter(n >= at_least_DAA)  %>%
    pull(Taxon)
  
  # Rank coefficients within each methodology
  DAA_ranks <- DAA_tibble %>% 
    filter(Taxon %in% which_taxa) %>% 
    group_by(DAA_tool) %>% 
    mutate(SignedRank = sign(coef) * rank(abs(coef), ties.method = 'first')) %>% 
    ungroup() %>% 
    select(Taxon, SignedRank, DAA_tool, meanRelAb)
  
  # Compute rank matrix
  wide_ranks <- DAA_ranks %>% 
    select(-meanRelAb) %>% 
    pivot_wider(names_from = c('DAA_tool'), 
                values_from = c('SignedRank')) 
  
  # Compute weights matrix
  weights_vector <- wide_ranks %>%
    left_join(distinct(DAA_ranks, Taxon, meanRelAb), by = "Taxon") %>%
    pull(meanRelAb)
  
  # Compute weighted correlation
  rank_matrix <- wide_ranks %>% select(where(is.numeric))
  n_tools <- ncol(rank_matrix)
  
  # outer applies a function to all pairs of indices from 1 to n_tools.
  # Thanks Gemini 2.5 pro for this snippet:
  weighted_cor_matrix <- outer(
    1:n_tools, # First vector for pairs (rows of the output matrix)
    1:n_tools, # Second vector for pairs (cols of the output matrix)
    
    # The function to apply to each pair of indices (i, j)
    # Vectorize() makes our function accept vectors of indices
    Vectorize(function(i, j) {
      weightedCorr(
        x = rank_matrix[[i]], 
        y = rank_matrix[[j]],
        method = "spearman",
        weights = weights_vector
      )
    })
  )
  dimnames(weighted_cor_matrix) <- list(colnames(rank_matrix), colnames(rank_matrix))
  return(weighted_cor_matrix)
}

# Extract mean relative abundances for dataset of interest
if(file.exists('Out/_Rdata/taxa_relAb_metrics.RDS')) {
  taxa_relAb_metrics.ls <- readRDS('Out/_Rdata/taxa_relAb_metrics.RDS')
} else {
  source('scripts/5.0_taxa_relab_metrics.R')
}

# Flatten the list in a large dataset
taxa_relab_metrics <- imap(taxa_relAb_metrics.ls, function(ds.ls, taxRank) {
  imap(ds.ls, function(dataset_tibble, dataset) {
    dataset_tibble %>% 
      mutate(taxRank = taxRank,
             dataset = dataset)
  }) %>% list_rbind()
}) %>% list_rbind()

# Subset input for a single dataset*CCE*taxrank
grouped_tibbles <- DAA %>%
  filter(taxRank == 'Genus'
         & DAA_tool %in% which_tools 
         & database %in% which_databases) %>%
  left_join(taxa_relab_metrics, by = c('Taxon', 'taxRank', 'dataset', 'database')) %>% 
  select(-taxRank) %>% 
  group_by(dataset, database) 

# Set names for each subgroup
group_names_df <- grouped_tibbles %>% 
  group_keys() %>% # to allow splitting by abundance type:
  left_join(CCE_metadata, 
            join_by('database' == 'Database')) %>% 
  unite(col = 'methodology', dataset, database, sep = '_')

tibble_list <- grouped_tibbles %>% 
  group_split(.keep = FALSE) 

tibble_list %<>% set_names(pull(group_names_df, methodology))

# This function splits the tibble list according to a method's characteristic
# and produces a separate heatmap. If no CCE_char is given, it produces a single
# heatmap for all.

split_tibList_by_samVar <- function(a_tibble_list, 
                                    at_least_DAA,
                                    CCE_char = FALSE) { # a column of CCE_metadata
  
  if(!c(CCE_char) %in% colnames(CCE_metadata) | CCE_char == FALSE) {
    factor_levels <- 'OneGroup'
    split_tibble_list <- list()
    split_tibble_list$OneGroup <- a_tibble_list
  } else { # Split the list by CCE_char
    group_names <- as.factor(group_names_df[[CCE_char]])
    split_tibble_list <- split(a_tibble_list, 
                               group_names)
    factor_levels <- levels(group_names)
  }
  
  for (category in factor_levels) {
    
    # Compute the spearman correlation across all CCE methods (and CCE subgroups if applicable)
    spearman_list <- map(
      split_tibble_list[[category]], 
      ~ rank_correlation_by_ds_db(.x, # a tibble
                                  at_least_DAA)) # 
    # Find the mean correlation
    avg_spearman <- Reduce('+', spearman_list) / length(spearman_list)
    
    # Reorder them
    #avg_spearman <- avg_spearman[which_tools,which_tools]
    
    # Plot
    heatmap_plot <- Heatmap(
      avg_spearman,
      col = colorRamp2(c(min(avg_spearman), (1+min(avg_spearman))/2.2, 1.0), 
                       c("#b45e2f", "#e4e8e1", "#00A759")), # Viridis palette,                
      show_heatmap_legend = FALSE, 
      row_names_gp = gpar(fontsize = 10),
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_column_names = FALSE,  
      show_row_dend = FALSE,
      row_names_side = "left",  
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          sprintf("%.2f", avg_spearman[i, j]), # Format number to 2 decimal places
          x, y, 
          gp = gpar(fontsize = 12) # Control font size here
        )
      },
      border = FALSE
      # A common alternative is to use a thin white border to separate cells:
      # rect_gp = gpar(col = "white", lwd = 1) 
    )
    
    pdf(
      file = paste0('Out/ISMB2025/clust_', category,'.pdf'), 
      width = 6, 
      height = 3.75,
      bg = 'white'
    )
    
    # 2. Draw the plot to the PDF file
    draw(heatmap_plot)
    
    # 3. CRITICAL: Close the device to save the file
    dev.off()
    
  } 
}

split_tibList_by_samVar(tibble_list, 7)
split_tibList_by_samVar(tibble_list, 'CCE_approach')
split_tibList_by_samVar(tibble_list, 'Taxonomy')

library(RColorBrewer) # For good color palettes

# --- Plot 1: The Average Correlation Heatmap (The Main Result) ---

# A diverging color palette is best for correlation
#color_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)



phylo_list <- lapply(cluster_list, as.phylo)
class(phylo_list) <- "multiPhylo"

consensus_tree <- consensus(phylo_list, p = 0.5, check.labels = TRUE)
plot(consensus_tree, main = "Majority-Rule Consensus of DA Method Clustering")
nodelabels(consensus_tree$node.label) # Add the support counts to the nodes

support_values <- as.numeric(consensus_tree$node.label) / length(phylo_list) * 100
plot(consensus_tree, main = "Majority-Rule Consensus (with % Support)")
nodelabels(paste0(round(support_values, 1), "%"))

for (i in seq(1,57)) {
  cluster <- cluster_by_ds_db(tibble_list[[i]], 5)
  plot(cluster, main = "Majority-Rule Consensus of DA Method Clustering")
  nodelabels(consensus_tree$node.label) # Add the support counts to the nodes
}


## Heatmap version 
