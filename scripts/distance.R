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
    dplyr::filter(Abundance.x*Abundance.y>0) %>% # subset to only ltaxa detected by both
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

##########################
## PCA on r.aitchison ####
##########################
tool_subset <- c('KB51', 'MOTUS', 'MPA_db2023', 'SM_gtdb_rs214_full', 'SM_genbank_202203')
# Create one dataframe by sample containing every tool community estimates
wide_abund.list <- melt_ps_list_glom(ps_family.ls, taxRank) %>% 
  dplyr::select(Sample, all_of(taxRank), Abundance, dataset, database) %>% 
  dplyr::filter(dataset == 'Feces' &
                  database %in% tool_subset) %>% 
  group_by(Sample) %>%
  pivot_wider(names_from = database, 
              values_from = Abundance, 
              values_fill = list(Abundance = 0)) %>%
  group_split()

# Compute pairwise distances between all tool pairs by sample
flat_dist <- map(wide_abund.list, compute_pairwise_dist)
stacked_matrix <- do.call(rbind, flat_dist)

# PCA on stacked matrix
pca_result <- prcomp(stacked_matrix, scale. = TRUE)  # Scaling to standardize the variance

# Extract PCA scores
pca_scores <- as.data.frame(pca_result$x)
prop_contrib <- summary(pca_result) %$% importance %>% .[2,1:2]*100 

#Loadings of each variable to the principal components
loadings <- as.data.frame(pca_result$rotation) %>% 
  rownames_to_column('tool_pair') %>% 
  select(PC1, PC2, tool_pair) %>% 
  tidyr::separate(tool_pair, into = c("tool1", "tool2"), sep = "_vs_")

# Visualize loadings for the first two components
ggplot(loadings, aes(x = PC1, y = PC2, colour = tool1, shape = tool2)) +
  geom_point(size=10) +
  #geom_text(hjust = 1.1, vjust = 1.1) +
  labs(title = "Principal component analysis of sample-wise distances between tool pairs",
       x = paste0('PC1 (',round(prop_contrib[1],0),'%)'),
       y = paste0('PC2 (',round(prop_contrib[2],0),'%)')) +
  theme_minimal()


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
                      na.value = 'white', name = 'distance') +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

# Violin plot
# dist_df %>% 
#   filter(dataset == 'Feces' &
#            tool1 %in% tool_subset & # avoid redundance
#            tool2 %in% tool_subset
#   ) %>%
#   pivot_longer(c('robustAitchison', 'bray'), names_to = 'dist') %>% 
#   ggplot(aes(x = tool2, y = value, fill = tool2)) +
#   geom_violin() +
#   facet_grid(dist ~ tool1, scales = 'free') +
#   theme(axis.text.x = element_blank(),
#         axis.title.x = element_blank()) +
#   ggtitle(paste0('Distribution of between-sample pairwise dissimilarity (',taxRank,'-level)'))
# 
