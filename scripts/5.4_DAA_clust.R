library(pacman)
p_load(tidyverse, magrittr, ape, RColorBrewer, ComplexHeatmap)
source("scripts/myFunctions.R")

DAA <- Sys.glob('Out/DAA/*.tsv') %>% 
  map(read_tsv, show_col_types = FALSE) %>% 
  list_rbind()

num_DAA_tools <- DAA$DAA_tool %>% unique %>% length
which_tools <- c('DESeq2', 'edgeR', 'Aldex2','corncob','MaAsLin2',  'ZicoSeq', 'ANCOMBC2'#, 'radEMU'
                 )
which_databases <- c('KB45', 'KB90', 'SM_RefSeq_20250528', 'KB45_GTDB', 'KB90_GTDB', 'SM_gtdb-rs220-rep', 
                     'MPA_db2023', 'MOTUS')

# Function to compute spearman correlation between 
cluster_by_ds_db <- function(DAA_tibble, at_least_DAA) {
  
  # Keep taxa detected by most methods
  which_taxa <- DAA_tibble %>% 
    group_by(Taxon) %>% 
    summarise(n = n()) %>% 
    filter(n >= at_least_DAA)  %>%
    pull(Taxon); length(which_taxa)
  
  # Rank coefficients within each methodology
  DAA_ranks <- DAA_tibble %>% 
    filter(Taxon %in% which_taxa) %>% 
    group_by(DAA_tool) %>% 
    mutate(SignedRank = sign(coef) * rank(abs(coef), ties.method = 'first')) %>% 
    ungroup() %>% 
    select(Taxon, SignedRank, DAA_tool)
  
  # Compute manhattan distance
  wide_ranks <- DAA_ranks %>% 
    pivot_wider(names_from = c('DAA_tool'), 
                values_from = 'SignedRank') #%>% 
  # select(where(is.numeric)) %>% t %>% 
  # dist(method = 'manhattan')
  
  wide_ranks %>% 
    select(where(is.numeric)) %>% 
    cor(method = 'spearman', use = "pairwise.complete.obs")
  
  #as.dist(1 - spearman_cor_mx)
  
  #dist_man_scaled <- dist_man/max(dist_man)
  #hclust(spearman_dist_mx, method = "ward.D2")
}

# Subset input for a single dataset*CCE*taxrank
grouped_tibbles <- DAA %>%
  filter(taxRank == 'Genus'
         & DAA_tool %in% which_tools 
         & database %in% which_databases) %>%
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
col_fun <- colorRamp2(c(0.3, 0.65, 1.0), c("#b45e2f", "#e4e8e1", "#00A759")) # Viridis palette

split_tibList_by_samVar <- function(a_tibble_list, 
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
    
    spearman_list <- map(split_tibble_list[[category]], 
                         ~cluster_by_ds_db(.x, 6))
    
    avg_spearman <- Reduce('+', spearman_list) / length(spearman_list)
    avg_spearman <- avg_spearman[which_tools,which_tools]
    heatmap_plot <- Heatmap(
      avg_spearman,
      
      col = col_fun,                
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
      
      # --- Borders ---
      # 'border_color = NA' in pheatmap removes the grid. The equivalent here is:
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

split_tibList_by_samVar(tibble_list)
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
