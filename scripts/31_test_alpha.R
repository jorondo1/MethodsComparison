library(pacman)
p_load(magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
       tidyverse, ggh4x, ggpubr, kableExtra,
       rstatix)

source("scripts/0_config.R")
source("scripts/30_test_plots.R")
theme_set(theme_light())
#
# $$$$$$$$\ $$$$$$\  $$$$$$\  $$\   $$\ $$$$$$$\  $$$$$$$$\       $$\   $$\
# $$  _____|\_$$  _|$$  __$$\ $$ |  $$ |$$  __$$\ $$  _____|      $$ |  $$ |
# $$ |        $$ |  $$ /  \__|$$ |  $$ |$$ |  $$ |$$ |            $$ |  $$ |
# $$$$$\      $$ |  $$ |$$$$\ $$ |  $$ |$$$$$$$  |$$$$$\          $$$$$$$$ |
# $$  __|     $$ |  $$ |\_$$ |$$ |  $$ |$$  __$$< $$  __|         \_____$$ |
# $$ |        $$ |  $$ |  $$ |$$ |  $$ |$$ |  $$ |$$ |                  $$ |
# $$ |      $$$$$$\ \$$$$$$  |\$$$$$$  |$$ |  $$ |$$$$$$$$\             $$ |
# \__|      \______| \______/  \______/ \__|  \__|\________|            \__|
#
# Alpha diversity tests ###

alpha_div <- read_rds('Out/_Rdata/alpha_div.RDS')

# Effect sizes ------------------------------------------------------------

# Rank-biserial correlation from the Mann-Whitney U (rstatix `statistic` = W of group A):
# 2U / (n1 n2) - 1 = P(A > B) - P(B > A)
wilcox.pdat <- alpha_div[['wilcox_tests']] %>%
  filter(Index %in% names(these_indices)) %>%
  mutate(rank_biserial = 2 * statistic / (n1 * n2) - 1) %>%
  method_rows_data() %>%
  split(.$Index)

# Plots -------------------------------------------------------------------

wilcox.plots <- imap(wilcox.pdat, function(df, idx) {
  method_rows_plot(
    df,
    x = rank_biserial,
    x_lab = paste0('Rank-biserial correlation (', these_indices[[idx]], ')'),
    zero_line = TRUE) +
    scale_x_continuous(
      breaks = scales::breaks_width(0.25),
      expand = expansion(add = 0.02))
}); wilcox.plots$Shannon

iwalk(wilcox.plots, function(plot, idx) {
  ggsave(paste0('Out/Manuscript/3.1.alpha_wilcox_', idx, '.pdf'),
         plot = plot,
         bg = 'white',
         width = 2000, height = 2000,
         units = 'px', dpi = 280)
})

# Tables ------------------------------------------------------------------

iwalk(wilcox.pdat, function(df, idx) {
  df %>%
    dplyr::select(Database, Dataset, rank_biserial, p) %>%
    mutate(
      value = paste0(format(rank_biserial, digits = 1),
                     " [p=", round(p, digits = 4), "]"),
      .keep = 'unused') %>%
    pivot_wider(names_from = 'Database', values_from = value) %>%
    kable() %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
    save_kable(paste0('Out/Manuscript/wilcox_', idx, '.html'))
})

# Example: violin plots ---------------------------------------------------

violin.pdat <- alpha_div[['plot_data']] %>%
  filter(Dataset == 'P19_Gut',
         Database %in% c('KB45', 'MPA_db2025'),
         Index == 'Shannon') %>%
  mutate(
    Group = factor(Grouping_var, labels = c('No', 'Yes')),  # A/B = diarr FALSE/TRUE
    Method = method_label(Database),
    Index = 'Shannon index')

violin.test <- violin.pdat %>%
  group_by(Index, Method) %>%
  wilcox_test(Index_value ~ Group) %>%
  add_xy_position(x = 'Group') %>%
  mutate(p_label = fmt_p(p))

violin.plot <- violin.pdat %>%
  ggplot(aes(x = Group, y = Index_value)) +
  geom_violin(
    aes(fill = Group),
    quantiles = 0.5,
    quantile.linetype = 1,
    linewidth = line_width) +
  geom_jitter(
    width = 0.15,
    size = point_size / 2,
    colour = 'grey20') +
  stat_pvalue_manual(
    violin.test,
    label = 'p_label',
    size = ex_text_size,
    bracket.size = line_width,
    tip.length = 0.02) +
  facet_wrap(~ Method) +
  scale_fill_manual(values = ex_fills, name = 'Diarrhea') +
  ex_theme

# Figure 4 ----------------------------------------------------------------

fig4 <- plot_grid(
  wilcox.plots$Shannon + tag_room,
  violin.plot + tag_room,
  ncol = 2,
  rel_widths = c(7, 2.3),
  labels = c('A.', 'B.'),
  label_size = 14,
  label_x = 0, hjust = 0); fig4

save_figure(fig4,
            pdf_path = 'Out/Manuscript/3.2.alpha_wilcox_example.pdf',
            eps_path = 'Out/Manuscript/PLOS/Fig4.eps')
