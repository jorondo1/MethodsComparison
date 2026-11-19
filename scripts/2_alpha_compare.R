pacman::p_load(magrittr, tidyverse, patchwork, ggh4x, ggridges, kableExtra, PairedData)
source('scripts/0_config.R')

#   $$\ $$\         $$$$$$$$\ $$$$$$\  $$$$$$\              $$\
#   $$ \$$ \        $$  _____|\_$$  _|$$  __$$\           $$$$ |
# $$$$$$$$$$\       $$ |        $$ |  $$ /  \__|          \_$$ |
# \_$$  $$   |      $$$$$\      $$ |  $$ |$$$$\             $$ |
# $$$$$$$$$$\       $$  __|     $$ |  $$ |\_$$ |            $$ |
# \_$$  $$  _|      $$ |        $$ |  $$ |  $$ |            $$ |
#   $$ |$$ |        $$ |      $$$$$$\ \$$$$$$  |$$\       $$$$$$\
#   \__|\__|        \__|      \______| \______/ \__|      \______|

# How much does each sample's alpha diversity move from one method to the next?
# Methods come in triads (alpha_triads, 0_config.R), each drawn as a closed loop (A -> B -> C -> A) so that
# every pair is adjacent once. Within a triad, only samples profiled by all 3 databases are kept.
# Differences are compared on indices scaled within dataset x database (median/MAD).

these_indices <- c(Richness = 'Richness', Shannon = 'Shannon', Hill_2 = 'Inverse Simpson', Tail = 'Tail')
axis_desc <- c(
  Richness = 'Richness',
  Shannon = 'Shannon entropy',
  Hill_2 = 'Inverse Simpson index',
  Tail = 'Tail statistic')

# One row per database: triad, position in the loop, labels
triads <- imap_dfr(alpha_triads, \(refs, type) imap_dfr(refs, \(dbs, ref)
  tibble(Type = type, Reference = ref, Database = names(dbs), Method = unname(dbs), pos = seq_along(dbs)))) %>%
  mutate(Triad = if_else(Type == 'Taxonomic abundances', 'Taxonomic', Reference),
         Panel = LETTERS[match(Triad, unique(Triad))])
triad_colours <- setNames(c('#2a78d6', '#eb6834', '#1baf7a'), unique(paste0(triads$Panel, '. ', triads$Triad)))

alpha_div <- read_rds('Out/_Rdata/alpha_div.RDS')[['plot_data']] %>%
  filter(Dataset %in% these_datasets, Database %in% triads$Database, Index %in% names(these_indices)) %>%
  group_by(Dataset, Database, Index) %>%
  # MAD is 0 when most samples share one value (Moss, sylph RefSeq: 64 % single-species): fall back to the SD
  mutate(z = (Index_value - median(Index_value)) / if_else(mad(Index_value) > 0, mad(Index_value), sd(Index_value)),
         pct = percent_rank(Index_value)) %>%
  ungroup()

# Samples profiled by all 3 databases of a triad
triad_dat <- alpha_div %>%
  inner_join(triads, by = 'Database') %>%
  group_by(Index, Triad, Dataset, Sample) %>% filter(n() == 3) %>%
  group_by(Index, Triad) %>% mutate(n = n_distinct(Dataset, Sample)) %>% ungroup()

# Loop segments 1 -> 2, 2 -> 3, 3 -> 1 (the last one ends on the repeated first axis, x = 4)
segments <- triad_dat %>%
  inner_join(triad_dat %>% select(Index, Triad, Dataset, Sample, pos_b = pos, Method_b = Method, value_b = Index_value, z_b = z),
             by = c('Index', 'Triad', 'Dataset', 'Sample'), relationship = 'many-to-many') %>%
  filter(pos_b == pos %% 3 + 1) %>%
  mutate(dz = abs(z_b - z),
         Strip = if_else(Type == 'Taxonomic abundances',
                         paste0(Panel, '. ', Type, ' (n = ', n, ')'),
                         paste0(Panel, '. ', Reference, ' (n = ', n, ')')))

# PANELS A-C : closed loops, one segment per sample and pair -------------
# Colour = dataset, alpha = |difference in scaled index|; strongest changes drawn last (on top)
loop_plot <- function(dat, idx, nested) {
  methods <- dat %>% distinct(pos, Method) %>% arrange(pos) %>% pull(Method)
  dat %>%
    arrange(dz) %>%
    ggplot(aes(x = pos, xend = pos + 1, y = Index_value, yend = value_b, colour = Dataset, alpha = dz)) +
    geom_segment(linewidth = 0.15) +
    { if (nested) facet_nested(~ Type + Strip) else facet_wrap(~ Strip) } +
    scale_x_continuous(breaks = 1:4, labels = c(methods, methods[1]), expand = expansion(add = 0.25)) +
    scale_colour_manual(values = dataset_colours, limits = these_datasets, labels = dataset_names) +
    scale_alpha_continuous(range = c(0.01, 0.8), limits = c(0, 3), oob = scales::squish,
                           breaks = 0:3, labels = c('0', '1', '2', '3+')) +
    labs(x = NULL, y = axis_desc[idx], colour = NULL, alpha = '|Difference in\nscaled index|') +
    guides(colour = guide_legend(nrow = 2, order = 1, override.aes = list(linewidth = 1, alpha = 1)),
           alpha = guide_legend(order = 2)) +
    theme_light(base_size = 8) +
    theme(!!!strip_theme,
          panel.spacing = unit(0.4, 'lines'),
          panel.grid.major.x = element_line(colour = 'grey70', linewidth = 0.3))
}

# PANEL D : largest pairwise difference per sample, by triad ------------
spread_plot <- function(dat) {
  spread <- dat %>%
    group_by(Panel, Triad, Dataset, Sample) %>%
    summarise(spread = max(dz), .groups = 'drop')
  x_max <- quantile(spread$spread, 0.99) # x axis cut at the 99th percentile (long tails, e.g. Inverse Simpson)
  spread %>%
    mutate(Triad = paste0(Panel, '. ', Triad),
           Host = if_else(Dataset %in% human_datasets, 'Human', 'Non-human'),
           Dataset = factor(Dataset, rev(these_datasets)),
           Strip = 'D. Per-sample disagreement') %>%
    ggplot(aes(x = spread, y = Dataset, fill = Triad)) +
    geom_density_ridges(
      stat = 'binline', binwidth = x_max / 30, boundary = 0, scale = 0.9, alpha = 0.5,
                        linewidth = 0.1, draw_baseline = FALSE) +
    facet_grid(Host ~ Strip, scales = 'free_y', space = 'free_y') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(labels = dataset_names, expand = expansion(add = c(0.1, 0.9))) +
    scale_fill_manual(values = triad_colours, guide = guide_legend(order = 3)) +
    coord_cartesian(xlim = c(0, x_max)) +
    labs(x = 'Largest |difference in scaled index| among the 3 pairs', y = NULL, fill = NULL) +
    theme_light(base_size = 8) +
    theme(!!!strip_theme, panel.spacing = unit(0.4, 'lines'))
}

alpha_figs <- imap(these_indices, function(label, idx) {
  dat <- filter(segments, Index == idx)
  loop_plot(filter(dat, Type == 'Sequence abundances'), idx, nested = TRUE) /
    (loop_plot(filter(dat, Type == 'Taxonomic abundances'), idx, nested = FALSE) + guides(colour = 'none') | # no Bee, Maple: legend from the top row
       spread_plot(dat)) +
    plot_layout(guides = 'collect') &
    theme(legend.position = 'bottom', legend.box = 'vertical', legend.margin = margin(0, 0, 0, 0))
})

# Figure 1 : Shannon
ggsave('Out/Manuscript/1.1.alpha_diff_Shannon.pdf', alpha_figs$Shannon,
       width = 190, height = 200, units = 'mm')
ggsave('Out/Manuscript/PLOS/Fig1.tiff', alpha_figs$Shannon, device = 'tiff', compression = 'lzw',
       bg = 'white', width = 190, height = 200, units = 'mm', dpi = 300)

# Figure SUPP : the same for the other indices
iwalk(alpha_figs[names(alpha_figs) != 'Shannon'], \(p, idx)
  ggsave(paste0('Out/Manuscript/1.2.alpha_diff_', idx, '.pdf'), p, width = 190, height = 200, units = 'mm'))

# Figure SUPP : pair matrix of the 9 databases ----------------------------
# Lower: within-dataset percentiles, one point per sample (diagonal = same rank with both databases).
# Upper: Spearman rho per dataset (bar = median). Samples profiled by both databases of each pair.
db_labels <- setNames(CCE_metadata$MethodName, CCE_metadata$Database)[triads$Database]

pair_matrix <- function(idx) {
  dat <- alpha_div %>% filter(Index == idx) %>% select(Dataset, Sample, Database, pct)
  pairs <- combn(triads$Database, 2, simplify = FALSE) %>%
    map_dfr(\(p) inner_join(filter(dat, Database == p[1]) %>% select(Dataset, Sample, x = pct),
                            filter(dat, Database == p[2]) %>% select(Dataset, Sample, y = pct),
                            by = c('Dataset', 'Sample')) %>%
              mutate(db_1 = factor(db_labels[p[1]], db_labels), db_2 = factor(db_labels[p[2]], db_labels)))
  rho <- pairs %>%
    group_by(db_1, db_2, Dataset) %>%
    summarise(rho = cor(x, y, method = 'spearman'), .groups = 'drop') %>%
    rename(db_1 = db_2, db_2 = db_1)   # rho goes above the diagonal
  ggplot() +
    geom_point(data = pairs, aes(x, y), size = 0.1, alpha = 0.2, colour = 'grey20') +
    geom_abline(data = distinct(pairs, db_1, db_2), aes(slope = 1, intercept = 0), colour = '#eb6834', linewidth = 0.3) +
    stat_summary(data = rho, aes(0.5, rho), fun = median, geom = 'crossbar', width = 0.9, linewidth = 0.25) +
    geom_point(data = rho, aes(0.5 + (as.numeric(factor(Dataset, these_datasets)) - 6.5) / 16, rho, colour = Dataset),
               size = 1.1) +
    facet_grid(db_2 ~ db_1, switch = 'y', drop = FALSE) +
    scale_colour_manual(values = dataset_colours, breaks = these_datasets, labels = dataset_names) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(min(0, rho$rho), 1)) +
    guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2.5))) +
    labs(x = NULL, y = NULL, colour = NULL,
         subtitle = paste0(these_indices[idx], '. Lower: percentile within dataset (diagonal = same rank). ',
                           'Upper: Spearman rho per dataset (bar = median).')) +
    theme_light(base_size = 7) +
    theme(!!!strip_theme, panel.spacing = unit(0.2, 'lines'), legend.position = 'bottom')
}

iwalk(these_indices, \(label, idx)
  ggsave(paste0('Out/Manuscript/1.3.alpha_pairs_', idx, '.pdf'), pair_matrix(idx), width = 250, height = 255, units = 'mm'))

# Pitman-Morgan tests -----------------------------------------------------
# Do two method pairs differ in the variance of their per-sample differences (scaled index)?
# Paired test on the samples common to both pairs; Holm across all indices and pairs of pairs.
# Sensitivity: only samples profiled by all 9 databases.

pair_diffs <- segments %>%
  transmute(Index, Dataset, Sample, diff = z - z_b,
            Pair = paste0(Triad, ': ', if_else(pos < pos_b, paste(Method, 'vs', Method_b), paste(Method_b, 'vs', Method))))

pitman_morgan <- function(diffs) {
  diffs %>%
    group_by(Index) %>%
    group_modify(\(d, k) {
      wide <- pivot_wider(d, id_cols = c(Dataset, Sample), names_from = Pair, values_from = diff)
      map_dfr(combn(unique(d$Pair), 2, simplify = FALSE), \(p) {
        sub <- wide %>% select(all_of(p)) %>% drop_na()
        tibble(pair_1 = p[1], pair_2 = p[2], n = nrow(sub),
               var_1 = var(sub[[1]]), var_2 = var(sub[[2]]),
               p = Var.test(sub[[1]], sub[[2]], paired = TRUE)$p.value)
      })
    }) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p, method = 'holm'))
}

all_9 <- alpha_div %>%
  group_by(Index, Dataset, Sample) %>% filter(n_distinct(Database) == nrow(triads)) %>%
  ungroup() %>% distinct(Index, Dataset, Sample)

pitman <- bind_rows(
  `All samples` = pitman_morgan(pair_diffs),
  `Samples in all 9 databases` = pitman_morgan(semi_join(pair_diffs, all_9, by = c('Index', 'Dataset', 'Sample'))),
  .id = 'Samples') %>%
  mutate(Index = these_indices[Index])

pitman %>%
  select(Samples, Index, pair_1, pair_2, n, var_1, var_2, p_adj) %>%
  mutate(across(c(var_1, var_2), ~ round(.x, 3)), p_adj = signif(p_adj, 3)) %>%
  kable() %>%
  kable_styling(bootstrap_options = c('striped', 'hover', 'condensed', 'responsive')) %>%
  save_kable('Out/Manuscript/alpha_pitman.html')

# p values: one pair of pairs per row, index x sample set as columns
pitman %>%
  mutate(column = paste0(Index, if_else(Samples == 'All samples', '', ' (all 9)'))) %>%
  pivot_wider(id_cols = c(pair_1, pair_2), names_from = column, values_from = p_adj) %>%
  kable(format = 'simple', digits = 3) %T>%
  print() %>%
  writeLines('Out/Manuscript/1.4.alpha_pitmanMorgan.txt')
