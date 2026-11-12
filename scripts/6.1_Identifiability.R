# Sample identifiability: share of each sample's signal a tool links to its reference database, on the tool's own
# scale (compare datasets and databases within a tool, not values between tools). Supplementary figure.
#   Kraken 2   reads Kraken assigns at species rank or below / all reads (Bracken summary, before redistribution)
#   sourmash   sum of f_unique_weighted (abundance-weighted share of the sample's k-mers in the gather genomes)
#   sylph      "% of reads detected in database" (genomes at > 95 % ANI)
#   MetaPhlAn  sum of the species rows (estimated reads of detected species / reads)
# mOTUs left out: its rate is a share of marker-gene reads only.

pacman::p_load(tidyverse, ggridges, scales)
source('scripts/0_config.R')
samples <- samples_qc()

ident <- bind_rows(
  read_rates('kraken_species', c('species', 'other', 'value')) %>%
    mutate(Panel = paste('Kraken', as.numeric(str_extract(Database, '\\d+')) / 100)),
  read_rates('sourmash', 'value') %>% mutate(Panel = 'sourmash'),
  read_rates('sylph', c('classified', 'unclassified', 'value')) %>% mutate(Panel = 'sylph'),
  read_rates('mpa', c('classified', 'unclassified', 'value')) %>% mutate(Panel = 'MetaPhlAn')
) %>%
  filter(Database %in% names(rate_databases), !is.na(value)) %>%
  semi_join(samples, by = c('Dataset', 'Sample')) %>%
  mutate(Reference = case_when(Panel == 'MetaPhlAn' ~ str_replace(Database, 'MPA_db', 'MetaPhlAn '),
                               str_detect(Database, 'GTDB|gtdb') ~ 'GTDB r220',
                               .default = 'RefSeq'))

panel_labels <- c(
  'Kraken 0.1'  = 'Kraken 2, conf. 0.10\nreads at species rank',
  'Kraken 0.45' = 'Kraken 2, conf. 0.45\nreads at species rank',
  'Kraken 0.9'  = 'Kraken 2, conf. 0.90\nreads at species rank',
  'sourmash'    = 'sourmash gather\nweighted k-mer share',
  'sylph'       = 'sylph\nreads detected (est.)',
  'MetaPhlAn'   = 'MetaPhlAn 4\nreads in species (est.)')
ref_colours <- c('RefSeq' = '#2a78d6', 'GTDB r220' = '#eb6834',
                 'MetaPhlAn 2022' = '#1baf7a', 'MetaPhlAn 2023' = '#eda100', 'MetaPhlAn 2025' = '#e87ba4')

ident %>%
  mutate(Panel = factor(Panel, names(panel_labels))) %>%
  ggplot(aes(x = value, y = Dataset, fill = Reference)) +
  geom_density_ridges(linewidth = 0.2, alpha = 0.6, scale = 1.2, rel_min_height = 0.005) +
  facet_grid(. ~ Panel, scales = 'free_x', labeller = as_labeller(panel_labels)) +
  scale_x_continuous(breaks = breaks_extended(3)) +
  scale_y_discrete(limits = rev(these_datasets), labels = dataset_names_n_qc(samples), expand = expansion(add = c(0.3, 1.3))) +
  scale_fill_manual(values = ref_colours, breaks = names(ref_colours)) +
  labs(x = "Share of the sample's signal linked to the database (tool-specific metrics)", y = NULL, fill = NULL) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_light(base_size = 9) +
  theme(strip.text = element_text(colour = 'grey20', size = 7),
        strip.background = element_rect(fill = 'grey92', colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        legend.position = 'bottom')

ggsave('Out/Manuscript/S.identifiability_ridges.pdf', width = 240, height = 110, units = 'mm')
