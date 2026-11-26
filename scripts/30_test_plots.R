# Helpers for the test figures (31_test_alpha.R, 32_test_beta.R):
# one row per dataset, one point per method; colour = configuration, shape = reference, size = p band.
# Source after scripts/0_config.R.

p_load(cowplot, colorspace)

# P-value bands -----------------------------------------------------------

p_levels <- c('n.s.', 'p < 0.05', 'p < 0.01')

p_band <- function(p) {
  cut(p,
      breaks = c(-Inf, 0.01, 0.05, Inf),
      labels = rev(p_levels),
      right = FALSE) %>%
    factor(levels = p_levels)
}

fmt_p <- function(p) {
  if_else(p < 0.001, 'p < 0.001', paste0('p = ', signif(p, 2)))
}

# Configurations ----------------------------------------------------------

# Tool + parameter, without the reference (shape shows it): "Kraken-Bracken 0.45", "sylph tax"
method_configs <- methods %>%
  arrange(Tool, Reference) %>%
  mutate(Config = case_when(
    Tool == 'sylph' ~ paste('sylph', if_else(Abundance == 'Sequence', 'seq', 'tax')),
    is.na(Param) ~ as.character(Tool),
    .default = paste(Tool, Param))) %>%
  dplyr::select(Database, Config)

# k shades of a tool colour, light to dark
shades <- function(col, k) {
  if (k == 1) return(col)
  amount <- seq(0.45, -0.35, length.out = k)  # > 0 lighten, < 0 darken
  map_chr(amount, \(a) if (a >= 0) lighten(col, a) else darken(col, -a))
}

config_colours <- methods %>%
  dplyr::select(Database, Tool) %>%
  inner_join(method_configs, by = 'Database') %>%
  arrange(Tool, Config) %>%
  distinct(Tool, Config) %>%
  mutate(colour = shades(tool_colours[[as.character(Tool[1])]], n()),
         .by = Tool) %>%
  { setNames(.$colour, .$Config) }

# Row layout --------------------------------------------------------------

ds_labels <- dataset_names_n()  # "Gut PD\n(n=719)", named by dataset
ds_order <- names(ds_labels)    # top to bottom

# Same offset for a method in every row; grouped by tool, then reference
offset_size <- 0.25
method_offset <- methods %>%
  arrange(Tool, Reference) %>%
  pull(Database) %>%
  { setNames(seq(offset_size, -offset_size, length.out = length(.)), .) }

ds_row <- function(ds) match(ds, rev(ds_order))

# Adds Tool, RefClass, Config, p band and y position to a Dataset x Database table with p
method_rows_data <- function(df) {
  df %>%
    filter(Dataset %in% ds_order,
           Database %in% methods$Database) %>%
    left_join(dplyr::select(methods, Database, Tool, RefClass),
              by = 'Database') %>%
    left_join(method_configs, by = 'Database') %>%
    mutate(
      Config = factor(Config, levels = names(config_colours)),
      p.signif = p_band(p),
      row = ds_row(Dataset),
      y = row + method_offset[Database])
}

# Main panel --------------------------------------------------------------

point_size <- 2
point_stroke <- 0.5
line_width <- 0.3
p_sizes <- setNames(c(1, 2, 3.2), p_levels)
refclass_solid <- c(`GTDB r220` = 17, RefSeq = 15, `Tool-specific` = 16)  # solid versions of refclass_shapes

# x: column of df (unquoted)
# shade_from: left edge of the row shading; use 0 on a log axis (log(-Inf) is NaN, log(0) = -Inf)
method_rows_plot <- function(df, x, x_lab, zero_line = FALSE, shade_from = -Inf) {
  shaded <- seq(1, length(ds_order), by = 2)

  ggplot(df, aes(x = {{ x }}, y = y)) +
    annotate(
      'rect',
      xmin = shade_from, xmax = Inf,
      ymin = shaded - 0.5, ymax = shaded + 0.5,
      fill = 'grey98') +
    { if (zero_line) geom_vline(
      xintercept = 0,
      colour = 'grey60',
      linewidth = line_width) } +
    geom_point(aes(colour = Config, shape = RefClass, size = p.signif)) +
    scale_colour_manual(
      values = config_colours,
      name = 'Method',
      guide = guide_legend(
        order = 1,
        nrow = 3, byrow = FALSE,  # one column per tool (sourmash shares sylph's)
        override.aes = list(shape = 16, size = 3))) +
    scale_shape_manual(
      values = refclass_solid,
      name = 'Reference database',
      guide = guide_legend(
        order = 2,
        override.aes = list(size = 3, colour = 'grey40'))) +
    scale_size_manual(
      values = p_sizes,
      drop = FALSE,
      name = 'p-value',
      guide = guide_legend(
        order = 3,
        override.aes = list(shape = 16, colour = 'grey40'))) +
    scale_y_continuous(
      breaks = seq_along(ds_order),
      labels = rev(unname(ds_labels)),
      expand = expansion(add = 0.2)) +  # outer offsets + a point radius
    labs(x = x_lab, y = NULL) +
    theme_light() +
    theme(
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = 'bottom',
      legend.box = 'vertical',
      legend.spacing.y = unit(8, 'pt'),  # between legends
      legend.key.height = unit(12, 'pt'),
      legend.key.spacing.y = unit(0, 'pt'),
      legend.margin = margin(0, 0, 0, 0))
}

# Example panel -----------------------------------------------------------

# Two groups, in factor-level order
ex_fills <- c('grey90', 'grey55')
ex_lines <- c('grey60', 'grey25')

ex_text_size <- 3  # mm, ~ axis text

ex_theme <- list(
  theme_light(),
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    legend.position = 'bottom',
    legend.title.position = 'top',
    !!!strip_theme))

# Figure ------------------------------------------------------------------

tag_room <- theme(plot.margin = margin(20, 5.5, 5.5, 5.5))  # space above panels for the A./B. labels

save_figure <- function(plot, pdf_path, eps_path) {
  ggsave(pdf_path, plot,
         bg = 'white',
         width = 2500, height = 2000,
         units = 'px', dpi = 260)
  ggsave(eps_path, plot,
         bg = 'white', device = 'eps',
         width = 2900, height = 2320,
         units = 'px')
}
