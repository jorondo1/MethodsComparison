library(pacman)
p_load(magrittr, tidyverse, RColorBrewer, plotly)

source('scripts/myFunctions.R')
source('scripts/5_DAA_fun.R')

compiled_DAA <- read_tsv('Out/DAA/Compiled_DAA_05.tsv')

# First a quick data manipulation function:
prep_data_for_heatmap <- function(df, count_by, which_values){
  # Fill 0 in nonexistent taxon/db combinations
  df %>% 
    pivot_wider(id_cols = where(is.character), 
                names_from = count_by, 
                values_from = all_of(which_values), 
                values_fill = 0) %>% 
    pivot_longer(cols = where(is.numeric), 
                 names_to = count_by, 
                 values_to = which_values) #%>% 
  #mutate(!!sym(which_values) := as.factor(!!sym(which_values)))
}

# 1. Look at how many DAA tools have identified a taxon 
# for each database, showing only taxa found by at least 2 DAA methods.
byatleast_DAA <- 3
sig_tax_DAA_count <- compiled_DAA %>% 
  filter(taxRank == 'Family' & dataset == 'NAFLD') %>%
  group_by(Taxon, database) %>% 
  summarise(count = n(), .groups = 'drop')  # count significant taxa per group

# keep only taxa detected by at least 2 methods
which_tax_DAA_count <- sig_tax_DAA_count %>% 
  dplyr::filter(count>=byatleast_DAA) %>% 
  pull(Taxon)

# Count DAA tools that identified each taxa in each dataset
count_tax_DAA <- sig_tax_DAA_count %>% 
  dplyr::filter(Taxon %in% which_tax_DAA_count) %>% 
  prep_data_for_heatmap(count_by = 'Taxon', which_values = 'count') %>% 
  mutate(count = as.factor(as.character(count))) 

# Plot !
count_tax_DAA %>% 
  ggplot(aes(x = database, y = Taxon)) +
  geom_tile(aes(fill = count), color = 'black', size =0.1) +
  scale_fill_manual(
    values = brewer.pal(n = length(levels(count_tax_DAA$count)), name="YlGnBu")
  ) + theme_minimal() +
  labs(fill = 'Number of\nDAA tools',
       title = 'Differential abundance agreement on taxon significance.',
       caption = 'Only taxa found to be differentially abundant by two or more approaches, at least once.') +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

ggsave('Out/DAA_count_DAA_NAFLD_F.pdf', bg = 'white', width = 2400, height = 1600, units = 'px', dpi = 180)

# 2. Same thing, but for each DAA_tool we want to know in how
# many database each taxon was found

sig_tax_db_count <- compiled_DAA %>% 
  filter(taxRank == 'Family' & dataset == 'NAFLD') %>%
  group_by(Taxon, DAA_tool) %>% 
  summarise(count = n(), .groups = 'drop')  # count significant taxa per group

# keep only taxa detected by at least 3 methods (filter out high edgeR+MPA false positive...)
which_tax_db_count <- sig_tax_db_count %>% 
  dplyr::filter(count>=byatleast_DAA) %>% 
  pull(Taxon)

# Count DAA tools that identified each taxa in each dataset
count_tax_db <- sig_tax_db_count %>% 
  dplyr::filter(Taxon %in% which_tax_db_count) %>% 
  prep_data_for_heatmap(count_by = 'Taxon', which_values = 'count') %>% 
  mutate(count = as.factor(as.character(count))) 

# Plot !
count_tax_db %>% 
  ggplot(aes(x = DAA_tool, y = Taxon)) +
  geom_tile(aes(fill = count), color = 'black', size =0.1) +
  scale_fill_manual(
    values = brewer.pal(n = length(levels(count_tax_db$count)), name="YlGnBu")
  ) + theme_minimal() +
  labs(fill = 'Number of composition estimation tools',
       title = 'Differential abundance agreement on taxon significance.',
       caption = 'Only taxa found to be differentially abundant by two or more approaches, at least once.') +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    legend.position = 'bottom',
    legend.title.position = 'top',
    legend.title = element_text(hjust = 0.5)
  ) + guides(fill = guide_legend(nrow = 1)) 

ggsave('Out/DAA_count_NAFLD_F.pdf', bg = 'white', width = 2400, height = 1600, units = 'px', dpi = 180)

# 3. Find a coefficient score (Cumulative effect sign index?)
# to represent the agreement in effect *sign*
sig_taxa_coef <- compiled_DAA %>% 
  filter(taxRank == 'Family' &
           dataset == 'NAFLD' & 
           !DAA_tool %in% c('ZicoSeq') ) %>% # ZicoSeq outputs all-positive coefficients...
  mutate(coef = case_when(coef > 0 ~ 1, 
                          coef < 0 ~ -1,
                          TRUE ~0)) %>% 
  group_by(Taxon, database) %>% 
  summarise(coef_score = sum(coef),
            .groups = 'drop') 

which_taxa_coef <- sig_tax_DAA_count %>% 
  dplyr::filter(count>1) %>% 
  pull(Taxon)

score_DAA <- sig_taxa_coef %>% 
  dplyr::filter(Taxon %in% which_taxa_coef) %>% 
  prep_data_for_heatmap(count_by = 'Taxon', which_values = 'coef_score') %>% 
  mutate(coef_score = factor(coef_score, levels = sort(unique(coef_score))))  

score_DAA %>% 
  ggplot(aes(x = database, y = Taxon)) +
  geom_tile(aes(fill = coef_score), color = 'black', size =0.1) +
  scale_fill_manual(
    values = brewer.pal(n = length(levels(score_DAA$coef_score)), name="PiYG")
  ) +
  # scale_fill_gradient2(low = 'darkgreen',
  #                      mid = 'white', 
  #                      high = 'darkred', 
  #                      midpoint = 0) +
  labs(fill = 'Cumulative Effect Sign Index',
       title = 'Frequency of significance calling for taxa found significant by at least 2 DAA tools overall.') +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(), 
    legend.position = 'bottom',
    legend.title.position = 'top',
    legend.title = element_text(hjust = 0.5)
  ) +guides(fill = guide_legend(nrow = 1)) 

ggsave('Out/DAA_coef_NAFLD_F.pdf', bg = 'white', width = 2400, height = 1600, units = 'px', dpi = 180)

# 2. For each database/DAA combination, how many taxa were found as sig?
# Restricted to taxa found by at least 2 DAA Tools overall.
taxa_found_multiple <- compiled_DAA %>% 
  dplyr::filter(dataset == 'NAFLD' & taxRank == 'Family') %>% 
  group_by(database, Taxon) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  dplyr::filter(count>1) %>% 
  pull(Taxon)

count_DB <- compiled_DAA %>% 
  filter(taxRank == 'Family' & 
           Taxon %in% taxa_found_multiple &
           dataset == 'NAFLD'
  ) %>% 
  group_by(database, DAA_tool) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  prep_data_for_heatmap(count_by = 'database') %>% 
  mutate(count = as.numeric(as.character(count)))

count_DB %>% 
  ggplot(aes(x = DAA_tool, y = database)) +
  geom_tile(aes(fill = count), color = 'black', size =0.1) +
  # we want 0s to be white, then have a gradient from 1 to max(count)
  scale_fill_gradientn(
    colors = c("white", "red", "blue"),  # Custom gradient
    values = scales::rescale(c(0, 1, max(count_DB$count))),  # Ensure 0 maps to white
    limits = c(0, max(count_DB$count)),  # Define the range of the scale
    na.value = "white"  # Optional: white for NA values if needed
  ) + theme_minimal() +
  labs(title = 'Number of significant taxa found, for taxa overall found significant by \nat least two DAA tools.') +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )


compiled_DAA %>% 
  filter(!DAA_tool %in% c('edgeR', 'DESeq2')) %>% 
  ggplot(aes(x = coef, y = log10(adj.p), colour = DAA_tool)) + 
  geom_point(size = 3) +
  theme_minimal() 

full_DAA <- rbind(
  read_tsv('Out/DAA/Maaslin2.tsv'),
  read_tsv('Out/DAA/AncomBC2.tsv'),
  read_tsv('Out/DAA/edgeR.tsv'), # Too many taxa, needs dealing with !
  read_tsv('Out/DAA/DESEq2.tsv'),
  read_tsv('Out/DAA/radEmu.tsv'),
  read_tsv('Out/DAA/Aldex2.tsv'),
  read_tsv('Out/DAA/ZicoSeq.tsv')
)

full_DAA %>% # Scale coefficients 
  group_by(taxRank, dataset, database, DAA_tool) %>% 
  mutate(coef_scaled = scale(coef)) %>% 
  # Choose subset (at least taxRank and dataset)
  filter(#!DAA_tool %in% c('edgeR', 'DESeq2') &
           taxRank == 'Genus' &
           dataset == 'NAFLD') %>% 
  ggplot(aes(x = coef_scaled, y = log10(adj.p), colour = database)) +
  geom_point(size = 0.6) +
  facet_grid(DAA_tool~database, scales = 'free') +
  geom_hline(aes(yintercept = log10(0.05), linetype = "p = 0.05"), color = "red", size = 0.5) +
  geom_hline(aes(yintercept = log10(0.01), linetype = "p = 0.01"), color = "blue", size = 0.5) +
  scale_linetype_manual(name = "", values = c("p = 0.05" = "dotted", "p = 0.01" = "dotted")) +
  scale_colour_manual(values = tool_colours) +
  guides(
    linetype = guide_legend(order = 1, override.aes = list(color = c("blue", "red"))),
    colour = guide_legend(order = 2),  # Keep the colour legend
    size = guide_legend(order = 4),    # Keep the size legend
    shape = guide_legend(order = 3)    # Keep the shape legend
) + theme(
  axis.text.y = element_blank(),
  axis.text.x = element_text(size = 5),
  axis.ticks = element_blank()
)

ggsave('Out/DAA_full_NAFLD_G.pdf', bg = 'white', width = 2400, height = 1600, units = 'px', dpi = 180)


fig <- plot_ly(
  data = count_tax_db,
  x = ~DAA_tool,
  y = ~Taxon,
  z = ~database,
  color = ~count,
  colors = colorRamp(c("blue", "green", "red")),
  type = "mesh3d"
)








