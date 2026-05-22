library(pacman)
p_load(magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
       tidyverse, kableExtra,
       igraph, ggraph, patchwork, pheatmap,
       rstatix)

source("scripts/0_Config.R")


#   $$\ $$\         $$$$$$$$\ $$$$$$\  $$$$$$\             $$$$$$\  
#   $$ \$$ \        $$  _____|\_$$  _|$$  __$$\           $$  __$$\ 
# $$$$$$$$$$\       $$ |        $$ |  $$ /  \__|          \__/  $$ |
# \_$$  $$   |      $$$$$\      $$ |  $$ |$$$$\            $$$$$$  |
# $$$$$$$$$$\       $$  __|     $$ |  $$ |\_$$ |          $$  ____/ 
# \_$$  $$  _|      $$ |        $$ |  $$ |  $$ |          $$ |      
#   $$ |$$ |        $$ |      $$$$$$\ \$$$$$$  |$$\       $$$$$$$$\ 
#   \__|\__|        \__|      \______| \______/ \__|      \________|
#                         

  ##################################
# PCoA procrustes visualisation ###
####################################

procrustes_tests <- read_rds('Out/_Rdata/procrustes.RDS')
these_databases <- c('MPA_db2023', 'MPA_db2022', 'MPA_db2025', 
                     'MOTUS3','MOTUS4', 
                     'KB10','KB10_GTDB', 'KB45', 'KB45_GTDB', 'KB90', 'KB90_GTDB',
                     'SM_gtdb-rs220-rep', 'SM_RefSeq_20250528')

these_datasets <- c('Moss', 'NAFLD','PD', 'Bee','AD_Skin', 'RA_Gut',
                    'P19_Gut', 'P19_Saliva')

# Heatmap with hclust ------------------------------------------------------
# One heatmap for all datasets

# Inspect higher pvalues; we expect significant between all datasets
procrustes_subset <- procrustes_tests %>%
  filter(db1 %in% these_databases,
         db2 %in% these_databases,
         ds %in% these_datasets) 

procrustes_subset %>% 
  filter(pval>0.01)

# We only plot significant correlations p<0.01
# because super low correlations is indicative of
# the tool not being able to detect the community in a meaningful way
db_pairs <- procrustes_subset %>% 
  filter(pval<0.01) %>% 
  group_by(db1, db2) %>%
  summarise(
    cor_tendency = median(cor),
    cor_dispersion = mad(cor, constant = 1), # Should it scale *1.4826 ?
    n_datasets = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(cor_tendency))

# long format
plot_data <- db_pairs %>%
  # unique database names
  bind_rows(
    db_pairs %>% 
      dplyr::select(
        db1 = db2, db2 = db1, 
        cor_tendency, cor_dispersion, n_datasets)
  ) %>%
  distinct() # no need for a diagonal, we want it to be NA

# Reorder factors in long df categories using hclust on values of a column;
# Built for 2-column factor (pairwise comparisons); db1 and db2 are pair ids

hclust_by_col <- function(data, colname){
  
  # All combinations
  all_dbs <- sort(unique(c(data$db1, data$db2)))
  
  # Build matrix
  cor_matrix <- matrix(
    1, 
    nrow = length(all_dbs), 
    ncol = length(all_dbs),
    dimnames = list(all_dbs, rev(all_dbs)))
  
  # Fill it
  for(i in 1:nrow(data)) {
    db1 <- as.character(data$db1[i])
    db2 <- as.character(data$db2[i])
    if(!is.na(pull(data,colname)[i])) {
      cor_matrix[db1, db2] <- pull(data,colname)[i]
    }
  }
  
  # hierarchical clustering
  hclust_result <- hclust(as.dist(1 - cor_matrix))
  ordered_dbs <- all_dbs[hclust_result$order]
  
  # reorder factors
  data %<>%
    mutate(
      db1 = factor(db1, levels = rev(ordered_dbs)),
      db2 = factor(db2, levels = ordered_dbs)
    ) 
  return(data)
}

# -- PLOT DATA Tendency (median) ---
cor_hclust.pdat <- plot_data %>% 
  hclust_by_col("cor_tendency") 

# Recode names
recode_names <- setNames(CCE_metadata$MethodName, CCE_metadata$Database)

cor_hclust_names.pdat <- cor_hclust.pdat %>% 
  mutate(db1 = recode(db1, !!!recode_names),
         db2 = recode(db2, !!!recode_names)) 

# Theme
theme_cor_mx <- list(
  theme_minimal(),
  theme(
    axis.text.x = element_text(size = 9, angle = 90, hjust =1),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    #aspect.ratio = 1
  )#,
  #coord_equal()
)


disp_cor.plot <- cor_hclust_names.pdat %>%
  # Plot !
  ggplot(aes(x = db1, y = db2, fill = cor_tendency)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", cor_tendency),
        color = is.na(cor_tendency)), 
    size = 3) +
  scale_color_manual(
    values = c("TRUE" = "white", "FALSE" = "black"), 
    guide = "none")+
  scale_fill_gradient(
    high = "#57a2ff", 
    #  mid = "#a0c83e", 
    low = "#c75740",
    # midpoint = median(plot_data$cor_tendency),
    #limits = c(0.6, 1),
    name = "Median Procrustes correlation",
    na.value = "white"
  ) +
  theme_cor_mx +
  theme(
    legend.position = 'bottom',
    legend.title.position = 'top'
  ); disp_cor.plot

ggsave('Out/Manuscript/2.1.procrustes_tendency_heatmap.pdf',
       plot = disp_cor.plot,
       bg = 'white', width = 2000, height = 2200, 
       units = 'px', dpi = 260)


# --- PLOT Dispersion ---
cor_hclust_names.pdat %>% 
  ggplot(aes(x = db1, y = db2, fill = cor_dispersion)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", cor_dispersion),
        color = is.na(cor_dispersion)), 
    size = 3) +
  scale_color_manual(
    values = c("TRUE" = "white", "FALSE" = "black"), 
    guide = "none")+
  scale_fill_gradient(
    high = "#c75740", 
    #  mid = "#a0c83e", 
    low = "#57a2ff",
    # midpoint = median(plot_data$cor_dispersion),
    name = "Median absolute deviation in Procrustes correlation",
    na.value = "white"
  ) +
  theme_cor_mx +
  theme(
    legend.position = 'bottom',
    legend.title.position = 'top'
  )

ggsave('Out/Manuscript/2.2.procrustes_disp_heatmap.pdf',
       bg = 'white', width = 2000, height = 2200, 
       units = 'px', dpi = 260)


# Ordination of correlations ----------------------------------------------

# Convert tibble to list by dataset and name list elements
procruste.list <- procrustes_subset %>% 
  group_by(ds) %>% 
  # Set names dynamically:
  (\(d) setNames(
    group_split(d, .keep = FALSE), 
    group_keys(d) %>% pull(1)))()

# --- matrix building helper function ---
matrix_builder <- function(procrustes, db1, db2, colname) {
  all_dbs <- sort(unique(c(
    as.character(procrustes[[db1]]),
    as.character(procrustes[[db2]])
  )))
  
  n <- length(all_dbs)
  
  # build correlation matrix
  cor_matrix <- matrix(1, nrow = n, ncol = n,
                       dimnames = list(all_dbs, all_dbs))
  
  for (i in 1:nrow(procrustes)) {
    v1  <- as.character(procrustes[[db1]][i])
    v2  <- as.character(procrustes[[db2]][i])
    val <- procrustes[[colname]][i]
    if (!is.na(val)) {
      cor_matrix[v1, v2] <- val
      cor_matrix[v2, v1] <- val   # symmetric fill
    }
  }
  return(cor_matrix)
}

# testing
matrix_builder(procruste.list$AD_Skin, 'db1', 'db2', 'cor')

# --- PCoA wrapper ---
pcoa_from_cor <- function(cor_mat, ds_name) {
  
  dist_mat <- as.dist(1 - cor_mat) # correlation = distance
  pcoa <- cmdscale(dist_mat, k = 2, eig = TRUE)
  
  tibble(
    Database = rownames(cor_mat),
    PC1 = pcoa$points[, 1],
    PC2 = pcoa$points[, 2],
    ds = ds_name,
    var_PC1 = pcoa$eig[1] / sum(pcoa$eig[pcoa$eig > 0]) * 100,
    var_PC2 = pcoa$eig[2] / sum(pcoa$eig[pcoa$eig > 0]) * 100
  )
}

# Run over all datasets 
pcoa_results <- imap_dfr(procruste.list, function(df, ds_name) {
  mat <- matrix_builder(df, "db1", "db2", "cor")
  pcoa_from_cor(mat, ds_name)
}) %>% 
  # Set names & levels
  mutate(ds = factor(recode(ds, !!!dataset_names), levels = dataset_names))

# Pull and define palettes by ref db
all_refdbs <- CCE_metadata %>%
  filter(Database %in% these_databases) %>%
  pull(short_alpha_2) %>% unique() %>% sort()

refdb_colors <- setNames(
  c("#823D51", tool_colours[['MetaPhlAn4']], tool_colours[['mOTUs4']], "#515A82"),
  all_refdbs
)

# PCoA Examples -----------------------------------------------------------

pcoa <- read_rds("Out/_Rdata/pcoa_noVST.ls.RDS")

pcoa_plot <- function(data.list, dataset, database, legend_title){
  
  eig1 <- format(data.list[[dataset]][[database]]$eig[1], digits = 2)
  eig2 <- format(data.list[[dataset]][[database]]$eig[2], digits = 2)
  #label_name = 'Gametophyte compartment'
  #label_colours = c('saddlebrown', 'darkgreen')
  
  data.list[[dataset]][[database]]$metadata %>% 
    mutate(Database = CCE_names[database]) %>% 
    ggplot(aes(x = PCo1, y = PCo2, 
               colour = !!sym(grouping_variable[[dataset]]))) +
    stat_ellipse(aes(fill = !!sym(grouping_variable[[dataset]])), 
                 geom = 'polygon',
                 alpha = 0.05) +
    geom_point(colour = 'black',
               stroke = 0.3,
               aes(fill = !!sym(grouping_variable[[dataset]])),
               shape = 21) +
    #scale_colour_manual(values = label_colours) +
    #scale_fill_manual(values = label_colours) +
    facet_grid(~Database)+
    labs(x = paste0("PCo1 [", eig1,"%]"),
         y = paste0("PCo2 [", eig2,"%]"),
         fill = legend_title, colour = legend_title
    ) +
    theme(
      axis.title.x = element_text(vjust=0),
      !!!strip_theme
    )
}

p1 <- pcoa_plot(pcoa, 'RA_Gut', 'KB45_GTDB', 'RA status')
p2 <- pcoa_plot(pcoa, 'RA_Gut', 'KB45', 'RA status')

library(cowplot)
# 1. Extract the shared legend from B or C
shared_legend <- get_legend(p1 + theme(
  legend.title.position = 'top',
  legend.position = "bottom",
))

# 2. Build B and C without legends, stacked with equal height
right_col <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  ncol = 1,
  align = "v",
  axis = "lr"
)

# 3. Add the shared legend below the right column
right_col_with_legend <- plot_grid(
  right_col,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.15)  # tune legend height
)

# 4. Combine A with the right column
final <- plot_grid(
  disp_cor.plot ,
  right_col_with_legend,
  labels = c('A.', 'B.'),
  ncol = 2,
  rel_widths = c(7,4),
  align = "h",
  axis = "t"             # align tops, let bottoms float freely
)

final


ggsave(paste0('Out/Manuscript/2.4.pcoa_example.pdf'),
       bg = 'white', width = 2200, height = 1600,
       units = 'px', dpi = 220)



# Per-dataset PCoA plots ---------------------------------------
plot_list <- pcoa_results %>%
  group_by(ds) %>%
  group_split() %>%
  setNames(pcoa_results %>% pull(ds) %>% unique() %>% sort()) %>%
  imap(function(df, ds_name) {
    
    # Per-facet axis labels with % variance explained
    x_lab <- sprintf("Eigenvector 1 [%.1f%%]", df$var_PC1[1])
    y_lab <- sprintf("Eigenvector 2 [%.1f%%]", df$var_PC2[1])
    
    df %>%
      left_join(CCE_metadata, by = "Database") %>%
      # Force all levels present in every plot — required for legend merging
      mutate(short_alpha_2 = factor(short_alpha_2, levels = all_refdbs)) %>%
      ggplot(aes(PC1, PC2, colour = short_alpha_2, label = MethodName)) +
      ggrepel::geom_label_repel(
        size = 2.5, 
        label.padding = 0.1,
        box.padding = 0.1,
        max.overlaps = 20,
        linewidth = 0,
        fill = NA,        # suppress label background (avoids ggrepel registering its own fill scale)
        label.size = NA,  # suppress label border
        colour = "black"
      ) +
      geom_point(size = 2) +
      # Suppress legend: a single shared legend is added manually below.
      # guides = "collect" in patchwork was unreliable when guide objects differed
      # subtly across plots, producing duplicate legends.
      scale_colour_manual(values = refdb_colors, drop = FALSE, guide = "none") +
      labs(subtitle = ds_name, x = x_lab, y = y_lab) +
      theme_bw() + 
      theme(
        axis.title = element_text(size = 9),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
      )
  })

# Reorder datasets (flush unused ones)
plot_list %<>% .[dataset_names] %>% compact()

# --- Standalone legend ---
# Some datasets are missing a few tools and that messes with patchwork's handling
# of legend collection. Bild dummy plot to guarantee a single legend.
# Horizontal layout with nrow = 1 enforced before extraction.

legend_plot <- ggplot(
  data.frame(short_alpha_2 = factor(all_refdbs, levels = all_refdbs)),
  aes(x = 1, y = short_alpha_2, fill = short_alpha_2)) +
  geom_point(size = 3, shape = 21, colour = "black") +
  scale_fill_manual(values = refdb_colors, name = "Reference database") +
  guides(fill = guide_legend(nrow = 4)) +
  theme_void() +
  theme(legend.position = "right", legend.direction = "vertical")

# extract legend
plot_list$Legend <- cowplot::get_legend(legend_plot)

# --- Patchwork ---
# Parentheses required: without them plot_layout() attaches only to legend_row
(wrap_plots(plot_list, ncol = 3)) 

ggsave('Out/Manuscript/2.3.procrustes_ordination.pdf',
       bg = 'white', width = 2200, height = 2200, 
       units = 'px', dpi = 240)




# Correlations into interpretable groups  ---------------------------------

# Deduplicate procruste; needs some gymnastic because of floating point incosistencies;
# procruses cor were computed twice per pair (AxB and BxA) and cor are slightly different
# procrustes_dedup <- procrustes_subset %>% 
#   rowwise() %>% 
#   mutate(unique_pair = paste(
#     min(as.character(db1), as.character(db2)), 
#     max(as.character(db1), as.character(db2)))
#   ) %>% 
#   group_by(ds, unique_pair) %>% 
#   slice(1) %>% 
#   ungroup()
# 
# 
# # Median Correlation when db1 and db2 are the same tool, but different database
# 
# procrustes_subset %>% 
#   filter(
#     str_detect(db1, "KB") & str_detect(db2, "KB")
#   )
# 
# # join to procruste data
# procrustes_classified <- procrustes_dedup %>%
#   # Add metadata to tool 1
#   left_join(CCE_metadata, by = c("db1" = "Database")) %>%
#   rename(method1 = CCE_approach, family1 = tool_family, database1 = refdb, taxonomy1 = Taxonomy) %>%
#   dplyr::select(ds, db1, db2, cor, pval, unique_pair, method1, family1, database1, taxonomy1) %>%
#   # add metadata to tool 2
#   left_join(CCE_metadata, by = c("db2" = "Database")) %>%
#   rename(method2 = CCE_approach, family2 = tool_family, database2 = refdb, taxonomy2 = Taxonomy) %>%
#   dplyr::select(ds, db1, db2, cor, pval, unique_pair,
#                 method1, family1, database1, taxonomy1,
#                 method2, family2, database2, taxonomy2) %>%
#   # Classify comparison types
#   mutate(comparison_type = case_when(
#     family1 == family2 & database1 == database2 ~ "Within tool, across parameter or database version",
#     family1 == family2 & database1 != database2 ~ "Same tool, different database (DNA-to-DNA only)",
#     method1 == method2 & family1 != family2 & database1 == database2 ~ "Same database, different tool (DNA-to-DNA only)",
#     method1 != method2 ~ "Cross-approach",
#     method1 == method2 & family1 != family2 ~ "Different tool & database (within approach)",
#     TRUE ~ "other [UNEXPECTED]"
#   ))
# 
# procrustes_classified %>% 
#   distinct(comparison_type, family1, family2) %>% 
#   arrange(comparison_type)
# 
# procrustes_classified %>%
#   group_by(comparison_type) %>%
#   summarise(
#     median_cor = median(cor),
#     MAD_cor = mad(cor, constant = 1),
#     n_pairs = n()
#   ) %>% 
#   arrange(median_cor)
# 
# 
