# distance-based analyses of abundance vectors
# Focusing only on taxa in common

library(pacman)
p_load(magrittr, tidyverse, phyloseq, rlang, vegan, ggridges)

ps_genus.ls <- read_rds("Out/ps_rare_genus.ls.rds") 
ps_family.ls <- read_rds("Out/ps_rare_family.ls.rds") 

# Pour une paire d'outils donnée, pour chaque échantillon on veut savoir si 
# le vecteur d'abondances relatives des espèces identifiées en commun
# est similaire.

# full join by sample&taxa, then only keep taxa detected by both
# transform to relab sample+tool-wise
# then, by sample, compute distance between both relab vectors

compile_distances <- function(df, tool_pair, taxRank) {
  
  compute_dist <- function(df, dist) {
    df %>% dplyr::select(relAb.x, relAb.y) %>% 
      as.matrix %>% t %>% # transpose as matrix
      vegdist(method = dist)
  }

  df %>%
    taxa_tool_pairs(tool_pair, taxRank) %>% # filter ds for tool pair
    #dplyr::filter(Abundance.x*Abundance.y>0) %>% # subset to only ltaxa detected by both
    mutate(across(where(is.numeric), ~ replace_na(., 0))) %>% 
     group_by(Sample) %>% 
    # Compute relative abundances by Sample 
    # note : bray needs this, but raitchison is agnostic to different scales
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
taxRank <- 'Order'
dist_df <- melt_ps_list_glom(ps_family.ls, taxRank) %>% 
  dplyr::select(Sample, all_of(taxRank), Abundance, dataset, database) %>% 
  apply_ds_toolpairs(compile_distances, taxRank) 

dist_df %<>% 
  mutate(dataset = factor(dataset, levels = c('Saliva', 'Feces', 'Moss'))) 

tool_subset <- c('KB51', 'MOTUS', 'MPA_db2023', 'SM_gtdb_rs214_full')


# Ridge plot 
dist_df %>% 
  filter(dataset %in% c('Feces', 'Saliva') &
           tool1 %in% tool_subset & 
           tool2 %in% tool_subset) %>%
  pivot_longer(c('robustAitchison', 'bray'), names_to = 'dist') %>% 
  # Creating the ridge plot
  ggplot(aes(x = value, y = tool1, fill = tool2)) +
  geom_density_ridges(scale = 0.9, alpha = 0.5) +
  facet_grid(dataset ~ dist, scales = 'free_x') +
  theme_light() + 
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.35))) +
  labs(
    title = paste0('Distribution of Between-Sample Pairwise Dissimilarity (',taxRank,'-level)'),
    x = "Distance or dissimilarity",
    y = "Tool"
  ) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"))

ggsave('Out/distances.pdf', bg = 'white', 
       width = 1600, height = 1600, units = 'px', dpi = 180)


# Violin plot

dist_df %>% 
  filter(dataset == 'Feces' &
           tool1 %in% tool_subset & # avoid redundance
           tool2 %in% tool_subset
         ) %>%
  pivot_longer(c('robustAitchison', 'bray'), names_to = 'dist') %>% 
  ggplot(aes(x = tool2, y = value, fill = tool2)) +
  geom_violin() +
  facet_grid(dist ~ tool1, scales = 'free') +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle(paste0('Distribution of between-sample pairwise dissimilarity (',taxRank,'-level)'))

# Sort the tool names to create a consistent combination
dist_df_nonredundant <- dist_df %>%
  rowwise() %>%
  mutate(tool_combination = paste(sort(c(tool1, tool2)), collapse = ".")) %>%
  ungroup() %>%
  # Group by dataset, Sample, and the sorted tool combination
  group_by(dataset, Sample, tool_combination) %>%
  # Keep only the first occurrence of each unique combination
  slice(1) %>%
  ungroup() %>%
  select(-tool_combination)


# Heatmap by mean within-sample distance with sd added
dist_df %>% 
  filter(dataset == 'Saliva') %>%
  group_by(tool1, tool2) %>% 
  summarise(mean_dist = mean(robustAitchison),
            sd_dist = sd(robustAitchison), .groups = 'drop') %>% 
  ggplot(aes(tool1, tool2, fill = mean_dist)) +
  geom_tile(color = 'white') +
  geom_text(aes(label = sprintf("%.3f",sd_dist)), color = 'white') +
  scale_fill_gradient(low = 'navyblue', high = 'red', 
                      na.value = 'white', name = 'distance', limits = c(1,7)) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

