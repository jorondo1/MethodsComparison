library(pacman)
p_load(
  # syntax
  magrittr, tidyverse, purrr, furrr, foreach, doParallel, parallel,
  # metrics and stats
  phyloseq, DESeq2, vegan, rstatix,
  # plotting :
  ggridges, ggbeeswarm2, patchwork)

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

#############
## PCoA ######
###############

# Create long dataframe for plotting pcoas
# all sample_data variables will end up in the df, expect lots of NA

# 1. compute distances for every dataset,
#    store in same structure list with one sublist per distance
compute_pcoa_wrap <- function(ps) {
  lapply(ps, function(ds) {
    lapply(ds, function(db) {
      list(
        `Bray-Curtis` = compute_pcoa(db, 'bray'),
        `Robust Aitchison` = compute_pcoa(db, 'robust.aitchison')
        )
      })
    })
  }

pcoa_genus.ls <- compute_pcoa_wrap(ps_genus.ls)
pcoa_species.ls <- compute_pcoa_wrap(ps_species.ls)

# PCoA compilation function
# extracts the sample_data tibble, adds dataset, database and distance names.
compile_pcoa <- function(ps, ds, db, dist) {
  ps[[dist]]$metadata %>% # generate tibble for this iteration
    rownames_to_column('Sample') %>% 
    mutate(dataset = ds,
           database = db,
           distance = dist) %>% 
    mutate(across(where(is.character), as_factor)) %>% 
    tibble
}

# 2. Compile over all datasets and distances to create a single dataframe
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

p1 <- plot_ordination_dist(pcoa_samdata, 'Saliva', 'Bray-Curtis', 'lostSmell') 
p2 <- plot_ordination_dist(pcoa_samdata, 'Saliva', 'Robust Aitchison', 'lostSmell') +
  theme(strip.text.x = element_blank())

p1 / p2 + plot_layout(guides = 'collect') & 
  theme(legend.position = "bottom")

ggsave('Out/pcoa_saliva_lostSmell.pdf', bg = 'white', 
       width = 2800, height = 1600, units = 'px', dpi = 240)

##################
### perMANOVA #####
####################

# iterate over pcoa results to test a single-factor permanova with the 
# diss.mx as a multivariate response 

iterate_permanova <- function(pcoa.ls, ds, vars) {
  ds.ls <- pcoa.ls[[ds]] 
  map(names(ds.ls), function(db) { # iterate over databases
    db.ls <- ds.ls[[db]]
    map(names(db.ls), function(dist) { # iterate over distances
      dist.mx <- db.ls[[dist]]$dist.mx
      samData <- db.ls[[dist]]$metadata
      
      # permanova
      formula <- as.formula(paste("dist.mx ~", paste(vars, collapse = " + ")))
      res <- adonis2(formula = formula, 
              permutations = 1000,
              data = samData,
              parallel = 8)
      
      # parse r2 and p for each explanatory variable 
      lapply(seq_along(vars), function(i) {
        tibble(
          database = db,              
          dist = dist,
          variable = vars[i],         
          R2 = res$R2[i],   
          p = res$`Pr(>F)`[i]) 
        
          }) %>% list_rbind
      }) %>% list_rbind
    }) %>% list_rbind
  }

permanova_df.ls <- list()
permanova_df.ls[['Feces']] <- iterate_permanova(pcoa_species.ls, 'Feces', c('group', 'sex', 'diarr', 'vacc', 'age'))
permanova_df.ls[['Saliva']] <- iterate_permanova(pcoa_species.ls, 'Saliva', c('group', 'sex', 'lostSmell', 'vacc', 'age'))
permanova_df.ls[['Moss']] <- iterate_permanova(pcoa_species.ls, 'Moss', c('Compartment', 'Host', 'SoilpH', 'SoilTemp', 'Compartment*Host'))

p_value_lines <- function() {
  list(
  geom_hline(aes(yintercept = log10(0.05), linetype = "p = 0.05"), color = "red"),
  geom_hline(aes(yintercept = log10(0.01), linetype = "p = 0.01"), color = "blue"),
  scale_linetype_manual(name = "", values = c("p = 0.05" = "dashed", "p = 0.01" = "dashed")),
  guides(
    linetype = guide_legend(order = 1, override.aes = list(color = c("red", "blue"))),
    colour = guide_legend(order = 2),  # Keep the colour legend
    size = guide_legend(order = 4),    # Keep the size legend
    shape = guide_legend(order = 3)    # Keep the shape legend
  )
  )
}

permanova_df.ls[['Saliva']] %>% 
  ggplot(aes(x = R2, y = log10(p), colour = database, shape = variable)) +
  geom_point(size = 5) + 
  scale_colour_manual(values = tool_colours) +
  p_value_lines() + theme_light() +
  labs(y = expression("log"[10]~"p-value"), x = expression("R"^2)) 

set.seed(1); permanova_df.ls[['Saliva']] %>%
  ggplot(aes(x = dist, y = log10(p), colour = database)) +
  geom_beeswarm(aes(size = R2, shape = dist), stroke = 1.5, spacing = 4, 
                method = 'swarm') +  # Use geom_beeswarm to avoid complete overlaps
  facet_grid(. ~ variable, scales = 'free_x') +
  labs(y = expression("log"[10]~"p-value"), 
       colour = 'Tool/database', 
       size = expression("R"^2), 
       shape = 'Distance metric',
       x = '') +
  p_value_lines() + plot_theme() + theme_light() + 
  scale_shape_manual(values = c(4,5)) +
  scale_size_continuous(range = c(0.1,10))+
  theme(axis.text.x = element_blank())

ggsave('Out/permanova_saliva_species.pdf', bg = 'white', 
       width = 2800, height = 1600, units = 'px', dpi = 240)


