# Dataset summary: samples with a non-empty profile per database, and read depth after QC.

pacman::p_load(tidyverse)
source('scripts/0_config.R')
samples <- samples_qc() # every post-QC sample

# Heatmap: share of post-QC samples with a non-empty profile (any level), per database and dataset ---------
# Non-empty: Bracken species reads, sourmash k-mer share, sylph / MetaPhlAn classified %, mOTUs assigned > 0.
# NA rates count as empty (10 PD mOTUs 3 profiles: rate rerun pending, see TODO).
profiled <- bind_rows(
  read_rates('bracken', 'value'),
  read_rates('sourmash', 'value'),
  read_rates('sylph', c('value', 'unclassified', 'rate')),
  read_rates('mpa', c('value', 'unclassified', 'rate')),
  read_rates('motus', c('value', 'unassigned', 'rate'))
) %>%
  filter(Database %in% names(rate_databases)) %>%
  semi_join(samples, by = c('Dataset', 'Sample')) %>%
  group_by(Dataset, Database) %>%
  summarise(n = sum(value > 0, na.rm = TRUE), .groups = 'drop') %>%
  left_join(count(samples, Dataset, name = 'n_qc'), by = 'Dataset') %>%
  mutate(pct = 100 * n / n_qc)

ggplot(profiled, aes(x = Dataset, y = Database, fill = pct)) +
  geom_tile(colour = 'white', linewidth = 0.5) +
  geom_text(aes(label = n, colour = pct > 50), size = 2.3) +
  scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = 'grey20'), guide = 'none') +
  scale_fill_gradient(low = '#cde2fb', high = '#1c5cab', limits = c(0, 100),
                      name = '% of post-QC samples\nwith a non-empty profile') +
  scale_x_discrete(limits = these_datasets, labels = dataset_names_n_qc(samples)) +
  scale_y_discrete(limits = rev(names(rate_databases)), labels = rate_databases) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('Out/Manuscript/S.profiled_samples.pdf', width = 190, height = 120, units = 'mm')

# Dataset table (html): post-QC samples (diversity comparisons), samples with metadata (= a value for the grouping
# variable: the tests), and depth of all post-QC samples in millions of read pairs ---------------------------------
dataset_table <- samples %>%
  group_by(Dataset) %>%
  summarise(post_qc = n(),
            depth = sprintf('%.1f ±%.1f [%.1f - %.1f]', mean(pairs) / 1e6, sd(pairs) / 1e6, min(pairs) / 1e6, max(pairs) / 1e6)) %>%
  left_join(count(samples_grouped(), Dataset, name = 'with_metadata'), by = 'Dataset') %>%
  arrange(factor(Dataset, these_datasets)) %>%
  transmute(Dataset = dataset_names[Dataset], post_qc, with_metadata, depth)

dataset_table %>%
  add_row(Dataset = 'Total', post_qc = sum(.$post_qc), with_metadata = sum(.$with_metadata), depth = '') %>%
  knitr::kable('html', escape = FALSE,
               col.names = c('Dataset', 'Samples (post-QC)', 'Samples with metadata', 'Mean ± sd MBP<sup>1</sup> [range]')) %>%
  c('<p><sup>1</sup> Sequence counts based on post-QC paired reads</p>') %>%
  writeLines('Out/dataset_summary/dataset_table.html')
