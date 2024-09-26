# distance-based analyses of abundance vectors
# Focusing only on taxa in common

library(pacman)
p_load(magrittr, tidyverse, phyloseq, rlang, vegan)

ps_genus.ls <- read_rds("Out/ps_rare_genus.ls.rds") 
ps_family.ls <- read_rds("Out/ps_rare_family.ls.rds") 

# Pour une paire d'outils donnée, pour chaque échantillon on veut savoir si 
# le vecteur d'abondances relatives des espèces identifiées en commun
# est similaire.

# full join by sample&taxa, then only keep taxa detected by both
# transform to relab sample+tool-wise
# then, by sample, compute distance between both relab vectors
taxRank <- 'Family'

compile_distances <- function(df, tool_pair, taxRank) {
  
  compute_dist <- function(df, dist) {
    df %>% dplyr::select(relAb.x, relAb.y) %>% 
      as.matrix %>% t %>% # transpose as matrix
      vegdist(method = dist)
  }
  
  df %>%
    taxa_tool_pairs(tool_pair, taxRank) %>% # filter ds for tool pair
     dplyr::filter(Abundance.x*Abundance.y>0) %>% # subset taxa detected by both
     group_by(Sample) %>% 
    # Compute relative abundances by Sample
     mutate(relAb.x = Abundance.x / sum(Abundance.x),
               relAb.y = Abundance.y / sum(Abundance.y),
            .keep = 'unused') %>%
     group_modify(~ { # Use groupings to compute sample dist between tool pairs
        # Compute distance
        bray <- compute_dist(.x, 'bray')
        roba <- compute_dist(.x, 'robust.aitchison')
        sampleName <- unique(.x$Sample)
        
        # Return two symmetric rows (for heatmap)
        tibble(Sample = c(sampleName, sampleName),
               tool1 = c(tool_pair[1], tool_pair[2]),
               tool2 = c(tool_pair[2], tool_pair[1]),
               bray = c(bray, bray),
               robustAitchison = c(roba,roba)
              )
        }
     )
}

# Compute Bray distances
bray_df <- melt_ps_list_glom(ps_family.ls, taxRank) %>% 
  dplyr::select(Sample, all_of(taxRank), Abundance, dataset, database) %>% 
  apply_ds_toolpairs(compile_distances, taxRank) 

bray_df %<>% 
  mutate(dataset = factor(dataset, levels = c('Saliva', 'Feces', 'Moss'))) %>% 
  filter(dataset != 'Moss')

# Heatmap by mean within-sample distance with sd added
bray_df %>% 
  filter(dataset == 'Feces') %>%
  group_by(tool1, tool2) %>% 
  summarise(mean_dist = mean(robustAitchison),
            sd_dist = sd(robustAitchison), .groups = 'drop') %>% 
  ggplot(aes(tool1, tool2, fill = mean_dist)) +
  geom_tile(color = 'white') +
  geom_text(aes(label = sprintf("%.3f",sd_dist)), color = 'white') +
  scale_fill_gradient(low = 'navyblue', high = 'red', 
                      na.value = 'white', name = 'Robust\nAitchison\ndistance') +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )



