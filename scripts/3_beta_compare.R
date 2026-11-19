pacman::p_load(magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
               tidyverse, kableExtra,
               igraph, ggraph, patchwork, pheatmap,
               rstatix, cowplot, ggrepel)

source("scripts/0_config.R") # methods, method_label(), reference_colours, tool_shapes

# helper for heatmap olrdering and PCoA
pairs_to_matrix <- function(data, value, missing = 1) {
  data <- mutate(data, across(c(db1, db2), as.character))
  dbs <- sort(unique(c(data$db1, data$db2)))
  m <- matrix(NA_real_, length(dbs), length(dbs), dimnames = list(dbs, dbs))
  m[cbind(data$db1, data$db2)] <- data[[value]]
  mt <- t(m)
  m <- ifelse(is.na(m), mt, ifelse(is.na(mt), m, (m + mt) / 2)) # average both directions, or take the one present
  m[is.na(m)] <- missing
  diag(m) <- 1
  m
}

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
these_databases <- methods$Database # all 17 methods of the study (0_config.R)


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

# Method order: hierarchical clustering on 1 - median correlation
db_order <- pairs_to_matrix(db_pairs, 'cor_tendency') %>%
  { rownames(.)[hclust(as.dist(1 - .))$order] }

# Both directions of every pair (no diagonal: left empty), in clustering order
plot_data <- db_pairs %>%
  bind_rows(rename(db_pairs, db1 = db2, db2 = db1)) %>%
  distinct(db1, db2, .keep_all = TRUE) %>%
  mutate(db1 = factor(db1, levels = rev(db_order)),
         db2 = factor(db2, levels = db_order))

# Heatmap of one column of plot_data; low/high: colours of the scale ends
cor_heatmap <- function(data, value, name, low, high) {
  ggplot(data, aes(x = db1, y = db2, fill = .data[[value]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = if_else(is.na(.data[[value]]), '', sprintf("%.2f", .data[[value]]))), size = 3) +
    scale_fill_gradient(low = low, high = high, name = name, na.value = "white") +
    scale_x_discrete(labels = method_label) +
    scale_y_discrete(labels = method_label) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 9, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = 'bottom',
      legend.title.position = 'top'
    )
}

# --- PLOT Tendency (median) ---
disp_cor.plot <- cor_heatmap(plot_data, 'cor_tendency', "Median Procrustes correlation",
                             low = "#c75740", high = "#57a2ff"); disp_cor.plot

ggsave('Out/Manuscript/2.1.procrustes_tendency_heatmap.pdf',
       plot = disp_cor.plot,
       bg = 'white', width = 2000, height = 2200,
       units = 'px', dpi = 260)

# --- PLOT Dispersion ---
cor_heatmap(plot_data, 'cor_dispersion', "Median absolute deviation in Procrustes correlation",
            low = "#57a2ff", high = "#c75740")

ggsave('Out/Manuscript/2.2.procrustes_disp_heatmap.pdf',
       bg = 'white', width = 2000, height = 2200,
       units = 'px', dpi = 260)

# PCoA Examples -----------------------------------------------------------

pcoa <- read_rds("Out/_Rdata/pcoa_noVST.ls.RDS")

pcoa_plot <- function(data.list, dataset, database, legend_title){
  
  eig1 <- format(data.list[[dataset]][[database]]$eig[1], digits = 2)
  eig2 <- format(data.list[[dataset]][[database]]$eig[2], digits = 2)
  
  data.list[[dataset]][[database]]$metadata %>%
    mutate(Database = method_label(database)) %>%
    ggplot(aes(x = PCo1, y = PCo2,
               colour = !!sym(grouping_variable[[dataset]]))) +
    stat_ellipse(aes(fill = !!sym(grouping_variable[[dataset]])),
                 geom = 'polygon',
                 alpha = 0.05) +
    geom_point(colour = 'black',
               stroke = 0.3,
               aes(fill = !!sym(grouping_variable[[dataset]])),
               shape = 21) +
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

# ggsave('Out/Manuscript/PLOS/Fig2.eps',
#        bg = 'white', device = 'eps', # transparency warning ok !
#        width = 2800, height = 2000, units = 'px')



#   $$\ $$\         $$$$$$$$\ $$$$$$\  $$$$$$\             $$$$$$\
#   $$ \$$ \        $$  _____|\_$$  _|$$  __$$\           $$ ___$$\
# $$$$$$$$$$\       $$ |        $$ |  $$ /  \__|          \_/   $$ |
#   \_$$  $$   |      $$$$$\      $$ |  $$ |$$$$\             $$$$$ /
#   $$$$$$$$$$\       $$  __|     $$ |  $$ |\_$$ |            \___$$\
# \_$$  $$  _|      $$ |        $$ |  $$ |  $$ |          $$\   $$ |
#   $$ |$$ |        $$ |      $$$$$$\ \$$$$$$  |$$\       \$$$$$$  |
#   \__|\__|        \__|      \______| \______/ \__|       \______/

# Ordination of correlations ----------------------------------------------
# Per dataset, a PCoA of the methods on 1 - Procrustes correlation. One encoding 
# per methods table (0_config.R): 
# colour = reference database,
# shape = tool, 
# filled = sequence abundance / open = taxonomic abundance
# label = the remaining parameter (Kraken confidence, MetaPhlAn database year, mOTUs version).

connect_tools <- c('Kraken') # thin line through each tool's parameter series (e.g. confidence 0.10 -> 0.45 -> 0.90)

procrustes_pcoa <- procrustes_subset %>%
  group_split(ds) %>%
  map_dfr(\(d) {
    m <- pairs_to_matrix(d, 'cor')
    p <- cmdscale(as.dist(1 - m), k = 2, eig = TRUE)
    var <- p$eig[1:2] / sum(p$eig[p$eig > 0]) * 100
    tibble(ds = d$ds[1], Database = rownames(m), PC1 = p$points[, 1], PC2 = p$points[, 2],
           strip = sprintf("%s [%.0f%% · %.0f%%]", dataset_names[d$ds[1]], var[1], var[2]))
  }) %>%
  left_join(methods, by = "Database") %>%
  mutate(strip = factor(strip, levels = unique(strip[order(match(ds, names(dataset_names)))])),
         ds_n = dataset_names_n(sep = ' ')[ds],
         fill = if_else(Abundance == 'Sequence', reference_colours[as.character(Reference)], 'white'),
         label = coalesce(Param, ''))

procrustes_ordination.plot <- ggplot(procrustes_pcoa, aes(PC1, PC2)) +
  geom_path(data = \(d) d %>% filter(Tool %in% connect_tools) %>% arrange(Param),
            aes(group = interaction(Tool, Reference), colour = Reference), linewidth = 0.3, show.legend = FALSE) +
  geom_point(aes(colour = Reference, shape = Tool, fill = fill), size = 2, stroke = 0.7) +
  geom_text_repel(
    seed = 1,
    aes(label = label), 
    size = 2.3, 
    colour = 'grey25', 
    box.padding = 0.15, 
    point.padding = 0.1,
    min.segment.length = Inf, 
    max.overlaps = Inf) +
  facet_wrap(~ strip, ncol = 4, scales = 'free') +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  scale_colour_manual(values = reference_colours, name = "Reference database") +
  scale_shape_manual(values = tool_shapes, name = "Tool") +
  scale_fill_identity() +
  guides(colour = guide_legend(override.aes = list(shape = 21, fill = reference_colours)),
         shape = guide_legend(override.aes = list(fill = 'grey60', colour = 'grey30'))) +
  labs(x = "Eigenvector 1", y = "Eigenvector 2",
       caption = paste0("Filled: sequence abundance; open: taxonomic abundance. ",
                        "Labels: Kraken confidence, MetaPhlAn database year, mOTUs version.\n",
                        "Strips: variance explained by eigenvectors 1 ; 2.")) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'bottom',
    legend.box = 'vertical',
    legend.spacing.y = unit(0, 'pt'),     # gap between legends
    legend.margin = margin(0, 0, 0, 0),
    plot.caption = element_text(hjust = 0, colour = 'grey30'),
    !!!strip_theme
  ); procrustes_ordination.plot

ggsave('Out/Manuscript/2.3.procrustes_ordination.pdf', procrustes_ordination.plot,
       bg = 'white', width = 2200, height = 2000,
       units = 'px', dpi = 240)

# ggsave('Out/Manuscript/PLOS/Fig3.eps', procrustes_ordination.plot,
#        bg = 'white', device = 'eps',
#        width = 2500, height = 2500, units = 'px')

# --- ANCHOR poster
# This figure was plotted with the old unfaceted (patchwork) plot 
# # see https://github.com/jorondo1/MethodsComparison/blob/6fcb5e209c50b3cdfa8d128374bdabec3d462a13/scripts/3_beta_compare.R
# procrustes_ordination.plot +
#   facet_wrap(~ ds_n, ncol = 3, scales = 'free') +
#   theme_light() +
#   theme(
#     panel.grid = element_blank(),
#     axis.ticks = element_blank(),
#     axis.text = element_blank(),
#     axis.title = element_text(size = 10),
#     legend.position = 'bottom',
#     legend.box = 'vertical',
#     legend.spacing.y = unit(0, 'pt'),     # gap between the two legends
#     legend.margin = margin(0, 0, 0, 0),
#     !!!strip_theme_dark
#   )
# 
# ggsave('Out/Manuscript/2.3.procrustes_ordination_ANCHOR.pdf',
#        bg = 'white', width = 2000, height = 2000,
#        units = 'px', dpi = 210)

