library(pacman)
p_load(magrittr, tidyverse, phyloseq,
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
# ...

### Parsing results
compiled_DAA <- rbind(
  read_tsv('Out/DAA/Maaslin2.tsv'),
  read_tsv('Out/DAA/AncomBC2.tsv'),
    read_tsv('Out/DAA/edgeR.tsv'), # Too many taxa, needs dealing with !
  # read_tsv('Out/DAA/DESEq2.tsv), #TDB
  # read_tsv('Out/DAA/radEmu1.tsv),
  read_tsv('Out/DAA/Aldex2.tsv')
) %>% 
  mutate(Taxon = fct_relevel(Taxon, sort(levels(Taxon))))

# First a quick data manipulation function:
prep_data_for_heatmap <- function(df, count_by){
    # Fill 0 in nonexistent taxon/db combinations
    pivot_wider(df, names_from =count_by, values_from = count, values_fill = 0) %>% 
    pivot_longer(cols = where(is.numeric), names_to = count_by, values_to = 'count') %>% 
    mutate(count = as.factor(count))
}

# 1. Look at how many DAA tools have identified a taxon 
# for each database, showing only taxa found by at least 2 DAA methods.
count_DAA <- compiled_DAA %>% 
  filter(taxRank == 'Genus' & dataset == 'NAFLD') %>% 
  group_by(Taxon, database) %>% 
  # count number of significant taxa per group
  summarise(count = n(), .groups = 'drop') %>% 
  # keep only taxa detected by at least 2 methods
  dplyr::filter(count>1) %>% 
  prep_data_for_heatmap(count_by = 'Taxon')

color_palette_DAA <- scales::colour_ramp(c("red", "blue"))(
  seq(0, 1, length.out = length(levels(count_DAA$count)) -1)
  )

# PLOT
count_DAA %>% 
  ggplot(aes(x = database, y = Taxon)) +
  geom_tile(aes(fill = count), color = 'black', size =0.1) +
  scale_fill_manual(
    values = c('white', color_palette_DAA)
  ) + theme_minimal() +
  labs(fill = 'Number of\nDAA tools',
       title = 'Taxa found significant by at least 2 DAA tools.') +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

# 2. For each database/DAA combination, how many taxa were found as sig?
# Restricted to taxa found by at least 2 DAA Tools 
count_DB <- compiled_DAA %>% 
  filter(taxRank == 'Family' & dataset == 'NAFLD') %>% 
  group_by(database, Taxon) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  prep_data_for_heatmap(count_by = 'database') %>% 
  mutate(count = as.numeric(as.character(count)))

count_DB %>% View
  ggplot(aes(x = DAA_tool, y = database)) +
  geom_tile(aes(fill = count), color = 'black', size =0.1) +
# we want 0s to be white, then have a gradient from 1 to max(count)
  scale_fill_gradientn(
    colors = c("white", "red", "blue"),  # Custom gradient
    values = scales::rescale(c(0, 1, max(count_DB$count))),  # Ensure 0 maps to white
    limits = c(0, max(count_DB$count)),  # Define the range of the scale
    na.value = "white"  # Optional: white for NA values if needed
  ) + theme_minimal() +
  labs(fill = 'Number of sig taxa tools') +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )


















