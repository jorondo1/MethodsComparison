library(pacman)
p_load(magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
       tidyverse, ggh4x)

source("scripts/0_config.R")
source("scripts/30_test_plots.R")
theme_set(theme_light())

#  $$$$$$$$\ $$$$$$\  $$$$$$\  $$\   $$\ $$$$$$$\  $$$$$$$$\       $$$$$$$\  
#  $$  _____|\_$$  _|$$  __$$\ $$ |  $$ |$$  __$$\ $$  _____|      $$  ____| 
#  $$ |        $$ |  $$ /  \__|$$ |  $$ |$$ |  $$ |$$ |            $$ |      
#  $$$$$\      $$ |  $$ |$$$$\ $$ |  $$ |$$$$$$$  |$$$$$\          $$$$$$$\  
#  $$  __|     $$ |  $$ |\_$$ |$$ |  $$ |$$  __$$< $$  __|         \_____$$\ 
#  $$ |        $$ |  $$ |  $$ |$$ |  $$ |$$ |  $$ |$$ |            $$\   $$ |
#  $$ |      $$$$$$\ \$$$$$$  |\$$$$$$  |$$ |  $$ |$$$$$$$$\       \$$$$$$  |
#  \__|      \______| \______/  \______/ \__|  \__|\________|       \______/ 
#
# PERMANOVA tests ###

permanova.ds <- read_rds('Out/_Rdata/permanova.ds.RDS')

# Plot --------------------------------------------------------------------

permanova.pdat <- method_rows_data(permanova.ds)

permanova.plot <- method_rows_plot(
  permanova.pdat,
  x = R2,
  x_lab = expression('PERMANOVA R'^2*' (Bray-Curtis)'),
  shade_from = 0) +
  scale_x_log10(
    breaks = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2),
    labels = scales::label_number(drop0trailing = TRUE),
    expand = expansion(mult = 0.03)); permanova.plot

ggsave('Out/Manuscript/4.beta_permanova.pdf',
       plot = permanova.plot,
       bg = 'white',
       width = 2000, height = 2000,
       units = 'px', dpi = 280)

# Example: PCoA -----------------------------------------------------------

pcoa.ex <- read_rds('Out/_Rdata/pcoa_noVST.ls.RDS')[['RA_Gut']][c('KB45_GTDB', 'KB45')]

pcoa.strips <- method_label(names(pcoa.ex))

# % variance on PCo1, PCo2
pcoa.var <- map_chr(pcoa.ex, function(pcoa) {
  var <- pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]) * 100
  sprintf('PCo1 %.0f%%, PCo2 %.0f%%', var[1], var[2])
})

pcoa.pdat <- imap_dfr(pcoa.ex, function(pcoa, db) {
  tibble(
    Database = db,
    PCo1 = pcoa$coordinates[, 1],
    PCo2 = pcoa$coordinates[, 2],
    Group = pcoa$metadata$Group)
}) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Group = factor(Group, levels = c('N', 'RA'), labels = c('Control', 'RA')),
    Method = factor(pcoa.strips[Database], levels = pcoa.strips))

pcoa.labels <- permanova.ds %>%
  filter(Dataset == 'RA_Gut',
         Database %in% names(pcoa.ex)) %>%
  mutate(
    Method = factor(pcoa.strips[Database], levels = pcoa.strips),
    label = paste0('R² = ', round(R2, 3), ', ', fmt_p(p), '\n', pcoa.var[Database]))

pcoa.plot <- pcoa.pdat %>%
  ggplot(aes(x = PCo1, y = PCo2)) +
  stat_ellipse(
    aes(colour = Group),
    linewidth = line_width) +
  geom_point(
    aes(fill = Group),
    shape = 21,
    size = point_size,
    stroke = point_stroke,
    colour = 'grey20') +
  geom_text(
    data = pcoa.labels,
    aes(label = label),
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.2,
    size = ex_text_size,
    lineheight = 0.9) +
  facet_wrap(~ Method, ncol = 1, scales = 'free') +  # stacked: near-square ordinations
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +  # room for the label
  scale_fill_manual(values = ex_fills, name = 'RA status') +
  scale_colour_manual(values = ex_lines, name = 'RA status') +
  labs(x = 'Bray-Curtis PCoA') +
  ex_theme +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text())

# Figure 5 ----------------------------------------------------------------

fig5 <- plot_grid(
  permanova.plot + tag_room,
  pcoa.plot + tag_room,
  ncol = 2,
  rel_widths = c(7, 2.3),
  labels = c('A.', 'B.'),
  label_size = 14,
  label_x = 0, hjust = 0); fig5

save_figure(fig5,
            pdf_path = 'Out/Manuscript/4.beta_permanova_example.pdf',
            eps_path = 'Out/Manuscript/PLOS/Fig5.eps')
