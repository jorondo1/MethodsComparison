library(pacman)
p_load(tidyverse, magrittr)

DAA <- rbind(
  read_tsv('Out/DAA_NAFLD//Maaslin2.tsv'),
  read_tsv('Out/DAA_NAFLD/AncomBC2.tsv'),
  read_tsv('Out/DAA_NAFLD/edgeR.tsv'), # Too many taxa, needs dealing with !
  read_tsv('Out/DAA_NAFLD/DESEq2.tsv'),
  read_tsv('Out/DAA_NAFLD/corncob.tsv'),
  #read_tsv('Out/DAA_NAFLD/radEmu.tsv'),
  read_tsv('Out/DAA_NAFLD/Aldex2.tsv')
  #read_tsv('Out/DAA_NAFLD/ZicoSeq.tsv')
)

num_DAA_tools <- DAA$DAA_tool %>% unique %>% length

DAA_subset <- DAA %>% 
  filter(#!DAA_tool %in% c('edgeR', 'DESeq2') &
    taxRank == 'Genus' &
      dataset == 'NAFLD' &
      database == 'KB20') 

# Rank coefficients within each methodology

which_taxa <- DAA_subset %>% group_by(Taxon) %>% 
  summarise(n = n()) %>% 
  filter(n == num_DAA_tools) %$% Taxon; length(which_taxa)

DAA_ranks <- DAA_subset %>% 
  filter(Taxon %in% which_taxa) %>% 
  group_by(# taxRank, dataset, 
           database, DAA_tool) %>% 
  mutate(#methodology = paste0(database, '_', DAA_tool),
         Rank = rank(coef, ties.method = 'first'),
         AdjRank = Rank - (max(Rank) + 1) / 2) %>% # adjusted manhattan
  ungroup %>% 
  select(Taxon, AdjRank, DAA_tool)

dist_man <- DAA_ranks %>% 
  pivot_wider(names_from = 'DAA_tool', 
              values_from = 'AdjRank') %>% 
  select(where(is.numeric)) %>% t %>% 
  dist(method = 'manhattan')

hc <- hclust(dist_man, method = "ward.D2")

# Plot dendrogram
plot(hc, main = "Hierarchical Clustering of Methods (Manhattan on Adjusted Ranks)")

# Within-database distance, then clustering all together?