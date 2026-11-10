# Sample identifiability: the share of each sample's signal that a tool links to its reference database,
# on that tool's own scale. Supplementary figure. Replaces the violins of 8_Classification_rates.R.
#
# The metrics differ in kind, so compare datasets and databases within a tool, not values between tools:
#   Kraken 2    reads at species level after Bracken / read pairs (Bracken moves reads from higher ranks
#               down to species, so this is close to the any-rank classified share)
#   sourmash    sum of f_unique_weighted: abundance-weighted share of the sample's k-mers covered by the
#               minimum set of genomes gather picks. A lower bound (exact k-mers, matches < 50 kb dropped)
#   sylph       "% of reads detected in database": estimated from the coverage of genomes at > 95 % ANI
#   MetaPhlAn   sum of the species rows: estimated reads of detected species / reads (marker coverage based)
#   mOTUs       assigned / (assigned + unassigned) marker inserts. Assigned includes meta and ext mOTUs
#               (no reference genome). mOTUs under 3 marker genes count in neither, so shallow samples lose signal.
# Depth: reads after KneadData (pairs x 2 + orphans), from its logs.
#
# Input: Out/classification_rates/*.tsv, written by scripts/src/classification_rates.sh (HPC).
# Output: Out/Manuscript/S.identifiability_*.pdf and Out/identifiability/identifiability_summary.tsv

library(pacman)
p_load(tidyverse, ggridges, patchwork, ggh4x, scales)
source('scripts/0_config.R')

rate_dir <- 'Out/classification_rates'
out_dir <- 'Out/identifiability'
dir.create(out_dir, showWarnings = FALSE)

# Datasets in the study (Olive is not). Display order: human gut, other human, animal, plant and soil.
dataset_order <- c('P19_Gut', 'NAFLD', 'PD', 'RA_Gut', 'P19_Saliva', 'AD_Skin', 'Bee', 'Moss', 'Maple',
                   'Grapevine', 'Oilfield')
study_datasets <- read_tsv('scripts/src/datasets.tsv', comment = '#', col_names = c('Dataset', 'tsv'),
                           show_col_types = FALSE) %>% pull(Dataset)
stopifnot(setequal(dataset_order, study_datasets))
ds_names <- c(dataset_names, Grapevine = 'Grapevine', Oilfield = 'Oilfield Soil')

# Columns of the rate files (no header); path = ./data/<Dataset>/<Methodology>/...
read_rates <- function(tool, columns) {
  read_tsv(file.path(rate_dir, paste0(tool, '_classification_rate.tsv')), col_names = columns,
           col_types = cols(path = 'c', .default = 'd'), na = 'NA') %>%
    mutate(Dataset = str_split_i(path, '/', 3),
           Methodology = str_split_i(path, '/', 4)) %>%
    filter(Dataset %in% study_datasets)
}

#### Depth ####
reads <- read_rates('reads', c('path', 'pairs', 'reads')) %>%
  mutate(Sample = basename(path)) %>%
  select(Dataset, Sample, pairs, reads)

#### Per-tool values ####
kraken <- read_rates('bracken', c('path', 'species_reads')) %>%
  mutate(Sample = str_split_i(path, '/', 5)) %>%
  left_join(reads, by = c('Dataset', 'Sample')) %>%
  mutate(value = species_reads / pairs,
         Panel = paste('Kraken', format(as.numeric(str_extract(Methodology, '\\d+')) / 100, nsmall = 2)),
         Database = if_else(str_detect(Methodology, 'GTDB'), 'GTDB r220', 'RefSeq'))

sourmash <- read_rates('sourmash', c('path', 'value')) %>%
  filter(Methodology %in% c('SM_RefSeq_20250528', 'SM_gtdb-rs220-rep', 'SM_gtdb-rs220-reps')) %>% # same database, folder name differs
  mutate(Sample = str_remove(basename(path), '_(RefSeq_\\d+|gtdb-reps?-rs220|gtdb-rs220-reps?)_gather\\.csv$'),
         Panel = 'sourmash',
         Database = if_else(str_detect(Methodology, 'RefSeq'), 'RefSeq', 'GTDB r220'))

sylph <- read_rates('sylph', c('path', 'classified', 'unclassified', 'value')) %>%
  filter(Methodology %in% c('SYLPH_RefSeq', 'SYLPH_gtdb_rs220')) %>% # glob_r232 and fungi_refseq: Maple only
  mutate(Sample = str_remove(basename(path), '_profile\\.log$'),
         Panel = 'sylph',
         Database = if_else(Methodology == 'SYLPH_RefSeq', 'RefSeq', 'GTDB r220'))

metaphlan <- read_rates('mpa', c('path', 'classified', 'unclassified', 'value')) %>%
  mutate(Sample = str_split_i(path, '/', 5),
         Panel = 'MetaPhlAn',
         Database = paste('MetaPhlAn', str_remove(Methodology, 'MPA_db')))

motus <- read_rates('motus', c('path', 'assigned', 'unassigned', 'value')) %>%
  mutate(Sample = str_remove(basename(path), '_profile\\.txt$'),
         Panel = 'mOTUs',
         Database = if_else(Methodology == 'MOTUS4', 'mOTUs 4', 'mOTUs 3'))

panel_order <- c('Kraken 0.10', 'Kraken 0.45', 'Kraken 0.90', 'sourmash', 'sylph', 'MetaPhlAn', 'mOTUs')
# Fixed categorical order (validated: dataviz palette, light mode, CVD-safe for adjacent pairs)
db_colours <- c('RefSeq' = '#2a78d6', 'GTDB r220' = '#eb6834',
                'MetaPhlAn 2022' = '#1baf7a', 'MetaPhlAn 2023' = '#eda100', 'MetaPhlAn 2025' = '#e87ba4',
                'mOTUs 3' = '#008300', 'mOTUs 4' = '#6250d6')

ident <- bind_rows(kraken, sourmash, sylph, metaphlan, motus) %>%
  filter(!is.na(value)) %>%
  select(Dataset, Sample, Panel, Database, value) %>%
  mutate(Panel = factor(Panel, panel_order),
         Database = factor(Database, names(db_colours)),
         Dataset = factor(Dataset, rev(dataset_order)))
stopifnot(!anyNA(ident$Panel), !anyNA(ident$Database), !anyNA(ident$Sample))

# Kraken's denominator: read pairs from KneadData equal Kraken's own total; check where the kreports remain
kraken_check <- read_rates('kraken', c('path', 'classified', 'unclassified', 'rate')) %>%
  mutate(Sample = str_split_i(path, '/', 5)) %>%
  inner_join(reads, by = c('Dataset', 'Sample')) %>%
  summarise(n = n(), mismatched = sum(classified + unclassified != pairs))
message('Kraken totals vs KneadData pairs: ', kraken_check$mismatched, ' of ', kraken_check$n, ' differ')

#### Summary table ####
summary_tbl <- ident %>%
  group_by(Panel, Database, Dataset) %>%
  summarise(n = n(), median = median(value), q25 = quantile(value, 0.25), q75 = quantile(value, 0.75),
            .groups = 'drop') %>%
  arrange(Panel, Database, desc(Dataset))
write_tsv(summary_tbl, file.path(out_dir, 'identifiability_summary.tsv'))

#### Figure A: identifiability per tool, with depth ####
# y labels: dataset name and sample count (KneadData)
n_lab <- reads %>% count(Dataset) %>%
  mutate(label = paste0(ds_names[Dataset], ' (', n, ')'))
y_labels <- setNames(n_lab$label, n_lab$Dataset)

metric_labels <- c(
  'Kraken 0.10' = 'Kraken 2 + Bracken\nconfidence 0.10',
  'Kraken 0.45' = 'Kraken 2 + Bracken\nconfidence 0.45',
  'Kraken 0.90' = 'Kraken 2 + Bracken\nconfidence 0.90',
  'sourmash' = 'sourmash gather\nweighted k-mer share',
  'sylph' = 'sylph\nreads detected (est.)',
  'MetaPhlAn' = 'MetaPhlAn 4\nreads in species (est.)',
  'mOTUs' = 'mOTUs\nassigned insert share')

ridge_theme <- list(
  theme_light(base_size = 9),
  theme(strip.text = element_text(colour = 'grey20', size = 7, lineheight = 0.9),
        strip.background = element_rect(fill = 'grey92', colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.x = element_text(size = 6.5),
        legend.position = 'bottom'))

p_depth <- reads %>%
  mutate(Dataset = factor(Dataset, rev(dataset_order)), panel = 'Reads after QC') %>%
  ggplot(aes(x = reads, y = Dataset)) +
  geom_density_ridges(fill = 'grey75', colour = 'grey25', linewidth = 0.2, scale = 0.9,
                      jittered_points = TRUE, position = position_points_jitter(height = 0, yoffset = -0.05),
                      point_shape = '|', point_size = 1.2, point_alpha = 0.5, point_colour = 'grey20') +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale()), breaks = 10^(4:9)) +
  scale_y_discrete(labels = y_labels, expand = expansion(add = c(0.2, 1))) +
  facet_grid(. ~ panel) +
  labs(x = 'Reads (log scale)', y = NULL) +
  ridge_theme

p_ident <- ident %>%
  ggplot(aes(x = value, y = Dataset, fill = Database)) +
  geom_density_ridges(colour = 'grey25', linewidth = 0.2, alpha = 0.6, scale = 0.9, rel_min_height = 0.005) +
  facet_grid(. ~ Panel, scales = 'free_x', labeller = as_labeller(metric_labels)) +
  scale_x_continuous(labels = label_number(accuracy = 0.1), breaks = breaks_extended(3)) +
  scale_y_discrete(labels = NULL, expand = expansion(add = c(0.2, 1))) +
  scale_fill_manual(values = db_colours, drop = TRUE) +
  labs(x = "Share of the sample's signal linked to the database (tool-specific scale, free axis per panel)",
       y = NULL, fill = 'Database') +
  guides(fill = guide_legend(nrow = 1)) +
  ridge_theme

fig_a <- p_depth + p_ident + plot_layout(widths = c(1.3, 7))
ggsave('Out/Manuscript/S.identifiability_ridges.pdf', fig_a, width = 280, height = 150, units = 'mm')
ggsave(file.path(out_dir, 'S.identifiability_ridges.png'), fig_a, width = 280, height = 150, units = 'mm',
       dpi = 200, bg = 'white')

#### Figure B: database effect within each tool, paired by sample ####
# Newer or larger database minus the other one, same sample and tool. Positive: the database links more.
pairs_spec <- tribble(
  ~Panel,        ~new,             ~old,             ~label,
  'Kraken 0.10', 'GTDB r220',      'RefSeq',         'Kraken 0.10\nGTDB - RefSeq',
  'Kraken 0.45', 'GTDB r220',      'RefSeq',         'Kraken 0.45\nGTDB - RefSeq',
  'Kraken 0.90', 'GTDB r220',      'RefSeq',         'Kraken 0.90\nGTDB - RefSeq',
  'sourmash',    'GTDB r220',      'RefSeq',         'sourmash\nGTDB - RefSeq',
  'sylph',       'GTDB r220',      'RefSeq',         'sylph\nGTDB - RefSeq',
  'MetaPhlAn',   'MetaPhlAn 2023', 'MetaPhlAn 2022', 'MetaPhlAn\n2023 - 2022',
  'MetaPhlAn',   'MetaPhlAn 2025', 'MetaPhlAn 2023', 'MetaPhlAn\n2025 - 2023',
  'mOTUs',       'mOTUs 4',        'mOTUs 3',        'mOTUs\n4 - 3 (also software)')

db_diff <- pairs_spec %>%
  pmap_dfr(function(Panel, new, old, label) {
    ident %>%
      filter(Panel == !!Panel, Database %in% c(new, old)) %>%
      mutate(Database = if_else(Database == new, 'new', 'old')) %>%
      pivot_wider(id_cols = c(Dataset, Sample), names_from = Database, values_from = value) %>%
      filter(!is.na(new), !is.na(old)) %>%
      mutate(diff = new - old, Comparison = label)
  }) %>%
  mutate(Comparison = factor(Comparison, pairs_spec$label))

write_tsv(db_diff %>% group_by(Comparison, Dataset) %>%
            summarise(n = n(), median_diff = median(diff), share_positive = mean(diff > 0), .groups = 'drop'),
          file.path(out_dir, 'identifiability_db_differences.tsv'))

fig_b <- db_diff %>%
  ggplot(aes(x = diff, y = Dataset)) +
  geom_vline(xintercept = 0, colour = 'grey50', linewidth = 0.3) +
  geom_density_ridges(fill = 'grey70', colour = 'grey25', linewidth = 0.2, scale = 0.9, rel_min_height = 0.005,
                      quantile_lines = TRUE, quantiles = 2) +
  facet_grid(. ~ Comparison, scales = 'free_x') +
  scale_y_discrete(labels = y_labels, expand = expansion(add = c(0.2, 1))) +
  scale_x_continuous(breaks = breaks_extended(3)) +
  labs(x = 'Paired difference in identifiability, same sample (vertical line in each ridge: median)', y = NULL) +
  ridge_theme

ggsave('Out/Manuscript/S.identifiability_db_differences.pdf', fig_b, width = 260, height = 140, units = 'mm')
ggsave(file.path(out_dir, 'S.identifiability_db_differences.png'), fig_b, width = 260, height = 140,
       units = 'mm', dpi = 200, bg = 'white')
