pacman::p_load(magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
               tidyverse, kableExtra,
               igraph, ggraph, patchwork, pheatmap,
               rstatix, cowplot, ggrepel)

source("scripts/0_config.R") # methods, method_label(), tool_colours, refclass_shapes

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

# DATASETS ----------------
procrustes_tests <- read_rds('Out/_Rdata/procrustes.RDS')
pcoa <- read_rds("Out/_Rdata/pcoa_noVST.ls.RDS")


#   $$\ $$\         $$$$$$$$\ $$$$$$\  $$$$$$\             $$$$$$\
#   $$ \$$ \        $$  _____|\_$$  _|$$  __$$\           $$  __$$\
# $$$$$$$$$$\       $$ |        $$ |  $$ /  \__|          \__/  $$ |
# \_$$  $$   |      $$$$$\      $$ |  $$ |$$$$\            $$$$$$  |
# $$$$$$$$$$\       $$  __|     $$ |  $$ |\_$$ |          $$  ____/
# \_$$  $$  _|      $$ |        $$ |  $$ |  $$ |          $$ |
#   $$ |$$ |        $$ |      $$$$$$\ \$$$$$$  |$$\       $$$$$$$$\
#   \__|\__|        \__|      \______| \______/ \__|      \________|
#
# PCoA procrustes visualisation ###

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

# Colour strips along the heatmap axes: 3-class database left of the y axis, approach under the x axis
strip_data <- methods %>%
  transmute(Database, RefClass, Approach = factor(approach_labels[Approach], levels = approach_labels))

# Lower-left triangle: median Procrustes correlation across datasets; upper-right: its median absolute deviation.
# legends: the legends to draw ('strips', 'gradients'), so that Fig 2 can lay them out on 2 rows (below)
cor_heatmap <- function(data, legends = c('strips', 'gradients')) {
  n <- nlevels(data$db1)
  data <- mutate(data, upper = as.integer(db1) + as.integer(db2) > n + 1) # diagonal: x + y = n + 1
  guide_if <- function(group, guide) if (group %in% legends) guide else 'none'
  triangle <- function(is_upper, value) {
    d <- filter(data, upper == is_upper)
    list(geom_tile(data = d, aes(fill = .data[[value]]), colour = 'white', linewidth = 0.5),
         geom_text(data = d, aes(label = sprintf('%.2f', .data[[value]])), size = 3))
  }
  ggplot(mapping = aes(x = db1, y = db2)) +
    triangle(FALSE, 'cor_tendency') +
    scale_fill_gradient(low = '#c75740', high = '#57a2ff', name = 'Median Procrustes correlation',
                        guide = guide_if('gradients', guide_colourbar(order = 3, theme = theme(legend.key.width = unit(2.8, 'cm'))))) +
    ggnewscale::new_scale_fill() +
    triangle(TRUE, 'cor_dispersion') +
    scale_fill_gradient(low = 'white', high = 'grey50', name = 'Median absolute deviation',
                        breaks = c(0, 0.05, 0.10), guide = guide_if('gradients', guide_colourbar(order = 4, theme = theme(legend.key.width = unit(2.8, 'cm'))))) +
    # strips sit in the panel just outside the tiles (rows and columns span 0.5 to n + 0.5)
    ggnewscale::new_scale_fill() +
    geom_tile(data = strip_data, aes(x = 0.1, y = Database, fill = RefClass), width = 0.4) +
    scale_fill_manual(values = refclass_colours, name = 'Database',
                      guide = guide_if('strips', guide_legend(order = 1))) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = strip_data, aes(x = Database, y = 0.1, fill = Approach), height = 0.4) +
    scale_fill_manual(values = approach_colours, name = 'Approach',
                      guide = guide_if('strips', guide_legend(order = 2))) +
    scale_x_discrete(labels = method_abbrev, expand = expansion(add = c(0, 0.1))) +
    scale_y_discrete(labels = method_abbrev, expand = expansion(add = c(0, 0.1))) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = 'bottom',
      legend.title.position = 'top',
      legend.title = element_text(hjust = 0)
    )
}

# Heatmap, then its legends in 2 columns: the 2 gradients, then the strips (Database, Approach), each stacked
heat_legend <- function(legends) {
  cowplot::get_legend(
    cor_heatmap(plot_data, legends) + 
      theme(
        legend.box = 'vertical',
        legend.box.just = 'left',
        legend.spacing.y = unit(4, 'pt')))
}
disp_cor.plot <- plot_grid(
  cor_heatmap(plot_data, legends = NULL),
  plot_grid(heat_legend('gradients'), heat_legend('strips'), nrow = 1),
  ncol = 1, rel_heights = c(1, 0.27)); disp_cor.plot

# PCoA Examples -----------------------------------------------------------

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
shared_legend <- cowplot::get_legend(p1 + theme(
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
  rel_widths = c(7,3),
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

connect_tools <- c('Kraken-Bracken') # thin line through each tool's parameter series (e.g. confidence 0.10 -> 0.45 -> 0.90)

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
         fill = if_else(Abundance == 'Sequence', tool_colours[as.character(Tool)], 'white'),
         label = coalesce(Param, ''))

procrustes_ordination.plot <- ggplot(procrustes_pcoa, aes(PC1, PC2)) +
  geom_path(data = \(d) d %>% filter(Tool %in% connect_tools) %>% arrange(Param),
            aes(group = interaction(Tool, Reference), colour = Tool), linewidth = 0.3, show.legend = FALSE) +
  geom_point(aes(colour = Tool, shape = RefClass, fill = fill), size = 2, stroke = 0.7) +
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
  scale_colour_manual(values = tool_colours, name = "Tool") +
  scale_shape_manual(values = refclass_shapes, name = "Reference database") +
  scale_fill_identity() +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 21, fill = tool_colours)),
         shape = guide_legend(order = 2, override.aes = list(fill = 'grey60', colour = 'grey30'))) +
  labs(x = "Eigenvector 1", y = "Eigenvector 2") +
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
