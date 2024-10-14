library(pacman)
p_load(magrittr, tidyverse, phyloseq, RColorBrewer,
       furrr, purrr, parallel,
       ANCOMBC, ALDEx2, radEmu, edgeR, DESeq2, Maaslin2)

source('scripts/myFunctions.R')
source('scripts/5_DAA_fun.R')

ps.ls <- read_rds('Out/ps_full.ls.rds')

out_path <- 'Out/DAA'
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

ncores <- detectCores() -1

# Aldex2
test_aldex <- compute_3_lvl(ps.ls, func = compute_aldex)
compile_3_lvl(test_aldex, func = compile_aldex) %>% 
  write_tsv('Out/DAA/Aldex2.tsv')

# ANCOM-BC2
test_ancombc2 <- compute_3_lvl(ps.ls, func = compute_ancombc2)
compile_3_lvl(test_ancombc2, func = compile_ancombc2) %>% 
  write_tsv('Out/DAA/AncomBC2.tsv')

# RadEmu
test_radEmu <- compute_3_lvl(ps.ls, func = compute_radEmu)
compile_3_lvl(test_radEmu, func = compile_radEmu) %>% 
  write_tsv('Out/DAA/radEmu1.tsv')

# MaAsLin2
capture_Maaslin_stdout <- compute_3_lvl(ps.ls, compute_Maaslin2)
compile_Maaslin(res_path = 'Out/DAA/Maaslin2/*/*/*/significant_results.tsv') %>% 
  write_tsv('Out/DAA/Maaslin2.tsv')

# EdgeR
test_edgeR <- compute_3_lvl(ps.ls, func = compute_edgeR)
compile_3_lvl(test_edgeR, func = compile_edgeR) %>% 
  write_tsv('Out/DAA/edgeR.tsv')

# DESeq2
test_DESeq2 <- compute_3_lvl(ps.ls, func = compute_DESeq2)
compile_3_lvl(test_DESeq2, func = compile_DESeq2) %>% 
  write_tsv('Out/DAA/DESeq2.tsv')

### Parsing results
compiled_DAA <- rbind(
  read_tsv('Out/DAA/Maaslin2.tsv'),
  read_tsv('Out/DAA/AncomBC2.tsv'),
  read_tsv('Out/DAA/edgeR.tsv'), # Too many taxa, needs dealing with !
  read_tsv('Out/DAA/DESEq2.tsv'),
  # read_tsv('Out/DAA/radEmu1.tsv),
  read_tsv('Out/DAA/Aldex2.tsv')
)

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
sig_taxa_count <- compiled_DAA %>% 
  filter(taxRank == 'Family' & dataset == 'NAFLD') %>%
  mutate(coef = case_when(coef > 0 ~ 1, coef < 0 ~ -1, TRUE ~0)) %>% 
  group_by(Taxon, database) %>% 
  # count number of significant taxa per group
  summarise(count = n(), coef_score = sum(coef),
            .groups = 'drop')

# keep only taxa detected by at least 2 methods
which_taxa <- sig_taxa_count %>% 
  dplyr::filter(count>1) %>% pull(Taxon)

# Count DAA tools that identified each taxa in each dataset
count_DAA <- sig_taxa_count %>% 
  dplyr::filter(Taxon %in% which_taxa) %>% 
  prep_data_for_heatmap(count_by = 'Taxon', which_values = 'count') %>% 
  mutate(count = as.factor(as.character(count))) 

color_palette_DAA <- scales::colour_ramp(c("white", "black"))(
  seq(0, 1, length.out = length(levels(count_DAA$count)))
  )

count_DAA %>% 
  ggplot(aes(x = database, y = Taxon)) +
  geom_tile(aes(fill = count), color = 'black', size =0.1) +
  scale_fill_manual(
    values = brewer.pal(n = length(levels(count_DAA$count)), name="YlOrRd")
  ) + theme_minimal() +
  labs(fill = 'Number of\nDAA tools',
       title = 'Taxa found significant by at least 2 DAA tools.') +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

ggsave('Out/DAA_count_NAFLD_F.pdf', bg = 'white', width = 2400, height = 1600, units = 'px', dpi = 180)

# PLOT with a coefficient score
score_DAA <- sig_taxa_count %>% 
  dplyr::filter(Taxon %in% which_taxa) %>% 
  prep_data_for_heatmap(count_by = 'Taxon', which_values = 'coef_score')


score_DAA %>% 
  ggplot(aes(x = database, y = Taxon)) +
  geom_tile(aes(fill = coef_score), color = 'black', size =0.1) +
  scale_fill_gradient2(low = 'blue',
                       mid = 'white', 
                       high = 'red', 
                       midpoint = 0) +
  labs(fill = 'Agreement score\nbetween DAA tools',
       title = 'Frequency of significance calling for taxa found significant by at least 2 DAA tools overall.') +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

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
  filter(DAA_tool != 'edgeR') %>% 
  ggplot(aes(x = coef, y = log10(adj.p), colour = database, shape = DAA_tool)) + 
  geom_point(size = 5)











