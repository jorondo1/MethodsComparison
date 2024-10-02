library(pacman)
p_load(magrittr, tidyverse, purrr, furrr, phyloseq, DESeq2, vegan, 
       rstatix, foreach, doParallel, patchwork)

Div_long <- read_rds('Out/Diversity_long.rds')
ps_species.ls <- read_rds("Out/ps_rare_species.ls.rds") 
ps_genus.ls <- read_rds("Out/ps_rare_genus.ls.rds") 

##############
## Alpha div ##
################

filter_and_add_samData <- function(df, ds, Rank, ps_list) {
  df %<>% filter(dataset == !!ds & Rank == !!Rank)
  
  # Extract sam_data from any phyloseq object of that dataset
  samData <- ps_list[[ds]][[1]]@sam_data %>% 
    as('data.frame') %>% 
    rownames_to_column('Sample')
  
  left_join(df, samData, by = 'Sample')
}

dataset <- filter_and_add_samData(
  Div_long, 'Saliva','Species',ps_species.ls
  ) %>% 
  group_by(database, index) 

dataset %>%
  shapiro_test(value) %>%
  filter(p<0.05)

# Check consistency of testing diversity difference between male/female
test_results <- dataset %>% 
  wilcox_test(value ~ lostSmell) %>% # conservative
  add_significance() %>% 
  select(database, index, p.signif) %>%
  mutate(p.signif = factor(p.signif, levels = c('ns','*','**','***', '****'))) # reorder factors

# Plot every database x index combination
test_results %>% # add the p-values to dataset
  left_join(dataset, join_by('database', 'index')) %>% 
  ggplot(aes(x = lostSmell, y = value, colour = p.signif)) +
  scale_color_discrete() +
  geom_boxplot(linewidth = 0.3) + 
  facet_grid(index~database, scales = 'free_y') +
  theme(
    panel.background = element_rect(fill = "grey95")  # Sets a lighter grey background
  )
ggsave('Out/alpha_saliva_lostSmell.pdf', bg = 'white', 
       width = 1800, height = 2400, units = 'px', dpi = 240)

###################
## Distance-based ##
#####################

# Create long dataframe for plotting pcoas
# all sample_data variables will end up in the df, expect lots of NA

# 1. compute distances for every dataset
pcoa_genus.ls <- lapply(ps_genus.ls, function(ds) {
  lapply(ds, function(db) {
    list(
      bray = compute_pcoa(db, 'bray'),
      robust.aitchison = compute_pcoa(db, 'robust.aitchison')
    )
  })
})

# PCoA compilation function
compile_pcoa <- function(ps, ds, db, dist) {
  ps[[dist]]$metadata %>% # generate tibble for this iteration
    rownames_to_column('Sample') %>% 
    mutate(dataset = ds,
           database = db,
           distance = dist) %>% 
    mutate(across(where(is.character), as_factor)) %>% 
    tibble
}

# 2. Compile!
pcoa_samdata <- iterate_distances(pcoa_genus.ls, compile_pcoa)

# Ordination plot between compartments :
plot_ordination_dist <- function(df, ds, dist, var) {
  df %>% 
    filter(dataset == ds
           & distance == dist
           & database %in% tool_subset) %>% 
    ggplot(aes(x = PCo1, y = PCo2, colour = !!sym(var))) + 
    stat_ellipse(level=0.9, geom = "polygon", alpha = 0.18, aes(fill = !!sym(var))) +   
    geom_point(size = 2) + 
    facet_grid(distance ~ database, scale = 'free')+
    theme(plot.title = element_text(size = 18),
          legend.title = element_text(colour="black", size=16, face="bold"),
          legend.text = element_text(colour="black", size = 14), 
          panel.spacing.x = unit(1, "lines"),
          axis.title = element_blank(),
          axis.text = element_blank(), 
          axis.ticks = element_blank()) + 
    guides(fill="none") 
}


p1 <- plot_ordination_dist(pcoa_full, 'Saliva', 'bray', 'lostSmell') 
p2 <- plot_ordination_dist(pcoa_full, 'Saliva', 'robust.aitchison', 'lostSmell') +
  theme(strip.text.x = element_blank())

p1 / p2 + plot_layout(guides = 'collect') & 
  theme(legend.position = "bottom")

ggsave('Out/pcoa_saliva_lostSmell.pdf', bg = 'white', 
       width = 2800, height = 1600, units = 'px', dpi = 240)

# simple perMANOVA comparison

