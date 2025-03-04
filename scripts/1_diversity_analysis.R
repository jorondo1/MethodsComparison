library(pacman)
p_load(
  # syntax
  magrittr, tidyverse, purrr, furrr, foreach, doParallel, parallel,
  # metrics and stats
  phyloseq, DESeq2, vegan, rstatix,
  # plotting :
  ggbeeswarm2, patchwork, grid, ggh4x)

source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/community_functions.R")
source('scripts/myFunctions.R')

#ps_species.ls <- read_rds("Out/ps_rare_species.ls.rds") 
#ps_genus.ls <- read_rds("Out/ps_rare_genus.ls.rds") 

ps.ls <- read_rds('Out/ps_full.ls.rds')
ps_rare.ls <- read_rds('Out/ps_rare.ls.rds')

#########################################
### Hill numbers and Tail diversity ####
#########################################

# Create long dataframe 
compile_diversity <- function(ps.ls) {
  
  # Compute diversity across indices
  div_rare.ls <- lapply(ps.ls, function(sublist) {
    lapply(sublist, div.fun, idx = c(0,1,2))
  })
  
  # Compile into 
  map(names(div_rare.ls), function(ds) { #iterate over dataset names
    ds_sublist <- div_rare.ls[[ds]] 
    map(names(ds_sublist), function(db) { # iterate over databases
      db_sublist <- ds_sublist[[db]]
      map(names(db_sublist), function(hill) { # iterate over index types
        values <- db_sublist[[hill]]
        tibble( # build dataset
          Sample = names(values),
          dataset = ds,
          database = db,
          index = hill,
          value = values
        ) %>% 
          mutate(across(where(is.character), as_factor))
      }) %>% list_rbind # collapse list into single df
    }) %>% list_rbind  
  }) %>% list_rbind
}

Div_long <- bind_rows(Species = compile_diversity(ps_rare.ls$Species), 
                      Genus = compile_diversity(ps_rare.ls$Genus), 
                      Family = compile_diversity(ps_rare.ls$Family),
                      .id = 'Rank') %>% 
  mutate(database = factor(database, levels = names(tool_colours)))

write_rds(Div_long, 'Out/Diversity_long_full.rds')

# Div_long <- read_rds('Out/Diversity_long.rds')
# Visualise differences in diversity across tools

Div_long %>% 
  filter(dataset=='NAFLD' & 
           index == 'H_1' &
           database %in% c('KB20', 'MOTUS','SM_genbank-2022.03', 'MPA_db2023') &
           Rank == 'Species') %>% 
  mutate(database = factor(database, levels = c('KB20', 'MOTUS','MPA_db2023', 'SM_genbank-2022.03'))) %>% 
  ggplot(aes(x = database, y = value, colour = database)) +
  geom_line(aes(group = Sample), alpha=0.3, linewidth = 0.1, colour = 'black') + 
  geom_point(size = 2) + 
  theme_light() + 
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  scale_colour_manual(values = tool_colours, labels = CCE_names)  +
  labs(y = 'Diversité', x = '', colour = "Méthode d'estimation\nde la composition")

ggsave('Out/diversity_filt_BFT708.pdf', bg = 'white', 
       width = 2300, height = 1300, units = 'px', dpi = 180)

Div_long %>% 
  filter(dataset!='AD_Skin' & 
    database != 'KB05' & !dataset %in% c('P19_Saliva', 'P19_Gut') &
           Rank != 'Family') %>% 
  # filter(Rank != 'Family' & index == 'H_0' & dataset != 'Feces') %>% 
  ggplot(aes(x = database, y = value, fill = database)) +
  geom_violin(size = 0.3) + 
  geom_line(aes(group = Sample), alpha=0.3, linewidth = 0.1) + 
  facet_nested(cols = vars(dataset, Rank), rows=vars(index), scales = 'free') +
  theme_light() + 
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values = tool_colours)

ggsave('Out/diversity_filt.pdf', bg = 'white', 
       width = 2400, height = 2400, units = 'px', dpi = 180)

##############
## Alpha div ##
################

plot_alpha_group <- function(ps.ls, ds, taxRank, Group) {
  dataset <- filter_and_add_samData(
    Div_long, ds, taxRank, ps.ls
    ) %>% 
    group_by(database, index) 

  # Check consistency of testing diversity difference between male/female
  test_results <- dataset %>% 
    wilcox_test(as.formula(paste("value ~", Group)), p.adjust.method = NULL) %>% # conservative
    add_significance() %>% 
    dplyr::select(database, index, p.signif) %>%
    mutate(p.signif = factor(p.signif, levels = c('ns','*','**','***', '****'))) # reorder factors
  
  # Plot every database x index combination
  p <- test_results %>% # add the p-values to dataset
    left_join(dataset, join_by('database', 'index')) %>% 
    ggplot(aes(x = !!sym(Group), y = value, colour = p.signif)) +
    scale_color_discrete() +
    geom_boxplot(linewidth = 0.3) + 
    facet_grid(index~database, scales = 'free_y') +
    theme(
      panel.background = element_rect(fill = "grey95")  # Sets a lighter grey background
    )
  
  ggsave(paste0('Out/alpha_',ds,'_',Group,'.pdf'), bg = 'white', plot = p,
       width = 1800, height = 2400, units = 'px', dpi = 240)
  
  return(p)
}
plot_alpha_group(ps_rare.ls, ds = 'P19_Saliva', 'Species', Group = 'lostSmell')

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

pcoa_genus.ls <- compute_pcoa_wrap(ps_rare.ls$Species)
pcoa_species.ls <- compute_pcoa_wrap(ps_rare.ls$Species)

# PCoA compilation function
# extracts the sample_data tibble, adds dataset, database and distance names.
compile_pcoa <- function(ps, ds, db, dist) {
  ps[[dist]]$metadata %>% # generate tibble for this iteration
    rownames_to_column('Sample') %>% 
    mutate(dataset = ds,
           database = db,
           distance = dist) %>% 
    mutate(across(where(is.factor), as.numeric)) %>% 
    tibble
}

# 2. Compile over all datasets and distances to create a single dataframe
pcoa_samdata <- iterate_distances(pcoa_species.ls, compile_pcoa)

# Ordination plot between compartments :
plot_ordination_dist <- function(df, ds, dist, var) {
  df %>% 
    dplyr::filter(dataset == ds
           & distance == dist
#           & database %in% tool_subset
           ) %>% 
    mutate(!!sym(var) := as.factor(!!sym(var))) %>% 
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

plot_beta_group <- function(pcoa_df, ds, Group) {
  p1 <- plot_ordination_dist(pcoa_df, ds, 'Bray-Curtis', Group) 
  p2 <- plot_ordination_dist(pcoa_df, ds, 'Robust Aitchison', Group) +
    theme(strip.text.x = element_blank())
  
  p <- p1 / p2 + plot_layout(guides = 'collect') & 
    theme(legend.position = "bottom")
  
  ggsave(paste0('Out/beta_',ds,'_',Group,'.pdf'), bg = 'white', plot = p,
         width = 3000, height = 1400, units = 'px', dpi = 240)
  return(p)
}
plot_beta_group(pcoa_samdata, 'NAFLD', 'Group')

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
              by = 'margin',
              na.action = na.exclude,
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
permanova_df.ls[['P19_Gut']] <- iterate_permanova(pcoa_species.ls, 'P19_Gut', c('sex', 'group', 'age'))
permanova_df.ls[['P19_Saliva']] <- iterate_permanova(pcoa_species.ls, 'P19_Saliva', c('sex', 'BMI', 'group', 'age'))
permanova_df.ls[['NAFLD']] <- iterate_permanova(pcoa_species.ls, 'NAFLD', c('Group', 'BMI'))
permanova_df.ls[['Moss']] <- iterate_permanova(pcoa_species.ls, 'Moss', c('Compartment', 'Host', 'SoilpH', 'SoilTemp'))
permanova_df.ls[['AD_Skin']] <- iterate_permanova(pcoa_species.ls, 'AD_Skin', c('Group', 'Gender', 'Ethnicity'))

write_rds(permanova_df.ls, 'Out/permanova_full.rds')

p_value_lines <- function() {
  list(
  geom_hline(aes(yintercept = log10(0.05), linetype = "p = 0.05"), color = "red"),
  geom_hline(aes(yintercept = log10(0.01), linetype = "p = 0.01"), color = "blue"),
  scale_linetype_manual(name = "", values = c("p = 0.05" = "dashed", "p = 0.01" = "dashed")),
  guides(
    linetype = guide_legend(order = 1, override.aes = list(color = c("blue", "red"))),
    colour = guide_legend(order = 2),  # Keep the colour legend
    size = guide_legend(order = 4),    # Keep the size legend
    shape = guide_legend(order = 3)    # Keep the shape legend
  )
  )
}

permanova_df.ls$P19_Gut %>% 
  ggplot(aes(x = R2, y = log10(p), colour = database, shape = variable)) +
  geom_point(size = 5) + 
  scale_colour_manual(values = tool_colours) +
  p_value_lines() + theme_light() +
  labs(y = expression("log"[10]~"p-value"), x = expression("R"^2)) 

permanova_df.ls$P19_Gut %>%
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

ggsave('Out/permanova_NAFLD_species.pdf', bg = 'white', 
       width = 2800, height = 1600, units = 'px', dpi = 240)
