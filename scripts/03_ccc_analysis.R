library(pacman)
p_load(magrittr, tidyverse, purrr, # syntax
       cccrm, cvequality, phyloseq, # specific
       patchwork, grid, ggh4x) # plotting

# Looking at the disagreement between community composition estimation tools
# when estimating various diversity indices from the resulting abundance
# matrices. This is a high-level exploratory analysis with the aim to 
# find DISagreement, as identical diversity does not prove identical compositions. 

source('scripts/myFunctions.R')

ps_raw.ls <- read_rds("Out/ps_raw.ls.rds")
ps_filt.ls <- read_rds("Out/ps_filt.ls.rds")
ps_species.ls <- read_rds("Out/ps_rare_species.ls.rds")
ps_genus.ls <- read_rds("Out/ps_rare_genus.ls.rds")
ps_family.ls <- read_rds("Out/ps_rare_family.ls.rds")

###################
### Sparseness ###
###################

# First, have a look at sparseness between filtered and unfiltered
# compute sparseness for all datasets
bind_rows(
  Unfiltered = compile_sparseness(ps_raw.ls), # id 1
  Filtered = compile_sparseness(ps_filt.ls), # id 2
  .id = 'filtered') %>% # flag id
  ggplot(aes(y = sparseness, x = database, fill = filtered)) +
  geom_col(position = 'dodge') + #theme_minimal() +
  facet_grid(dataset ~database, space = 'free', scales = 'free_x') +
  theme(axis.text.x = element_blank()) + ylim(c(0,1))

ggsave('Out/sparseness_filtering.pdf', bg = 'white', 
       width = 2400, height = 1600, units = 'px', dpi = 180)

###########################################
### Concordance Correlation Coefficient ####
###########################################

# Add tool category (DNA-to-DNA or DNA-to-Marker)
Div <- Div_long %>%
  mutate(across(where(is.factor), as.character)) %>% 
  mutate(type = case_when(
    str_starts(database, '^KB|^SM') ~ 'DNA',
    str_starts(database, '^MPA|^MOTUS') ~ 'marker'
  ))

# Apply the pairwise CCC calculation for each dataset and index group
ccc_pairwise_df <- Div %>%
  group_by(dataset, index, Rank) %>%
  group_modify(~ {
    # Create unique tool pairs within the group
    tools <- unique(.x$database)
    tool_pairs <- c(combn(tools, 2, simplify = FALSE), 
                    lapply(tools, function(x) c(x, x)))
    # For each instance of that pair, apply the cccvc_compile function
    map_dfr(tool_pairs, function(pair) cccvc_compile(.x, pair))
  }) %>%
  ungroup() %>% 
  mutate(CCC = case_when(tool1 == tool2 ~ NA,
                         TRUE ~ CCC) )

##################
##### HEATMAP ###
################

# Plot a heatmap for a given combination of dataset and index
plot_heatmaps <- function(df, dataset, index) {
  filtered <- df %>% 
    dplyr::filter(
      dataset == !!dataset, 
      index == !!index 
    ) 
  # Plot
  ggplot(filtered, aes(tool1, tool2, fill = CCC)) +
    geom_tile(color = 'white') +
    scale_fill_gradient(low = 'deeppink', high = 'navyblue', 
                        na.value = 'white', name = 'CCC',
                        limits = c(0,1)) +
    theme_minimal() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank()
    )
  }

# Make a heatmap grid for a given taxonomic rank by iterating over 
# each dataset & index combination 
heatmap_grid <- function(df, Rank) {
  ccc <- df %>% filter(Rank == !!Rank)
  
  combinations <- ccc %>% distinct(dataset,index)
  
  plots <- list()
  # Generate & compile plots
  for (i in seq_along(combinations$dataset)) {
    dataset <- combinations$dataset[i]
    index <- combinations$index[i]
    
    # Create the plot
    p <- plot_heatmaps(ccc, dataset, index)
    
    # Is plot first in row? 
    if (i %% 4 != 1) {
      # Remove the y-axis if not
      p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    } else {
        p <- p + guides(fill = "none")
    }
    plots[[i]] <- p
  }
  
  # Create row labels
  row_labels <- lapply(unique(combinations$dataset), function(label) {
    grid::textGrob(label, rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
  })
  
  # Create column labels
  col_labels <- lapply(unique(combinations$index), function(label) {
    grid::textGrob(label, gp = gpar(fontsize = 12, fontface = "bold"))
  })
  
  # Combine the plots with empty slots for labels
  wrap_plots(
    plot_spacer(), col_labels[[1]], col_labels[[2]], col_labels[[3]], col_labels[[4]],
    row_labels[[1]], plots[[1]], plots[[2]], plots[[3]], plots[[4]],
    row_labels[[2]], plots[[5]], plots[[6]], plots[[7]], plots[[8]],
    row_labels[[3]], plots[[9]], plots[[10]], plots[[11]], plots[[12]],
    row_labels[[4]], plots[[13]], plots[[14]], plots[[15]], plots[[16]],
    ncol = 5, nrow = 5,
    heights = c(0.2, 1, 1, 1,1), # Adjust height of the top row for column labels
    widths = c(0.2, 1, 1, 1, 1,1) # Adjust width of the first column for row labels
  ) + 
    plot_layout(guides = 'collect') & 
    theme(legend.position = "right")
}
heatmap_grid(ccc_pairwise_df, 'Species') 
ggsave('Out/ccc_Species.pdf', bg = 'white', width = 2400, height = 1600, units = 'px', dpi = 180)

heatmap_grid(ccc_pairwise_df, 'Genus')
ggsave('Out/ccc_Genus.pdf', bg = 'white', width = 2400, height = 1600, units = 'px', dpi = 180)

heatmap_grid(ccc_pairwise_df, 'Family')
ggsave('Out/ccc_Family.pdf', bg = 'white', width = 2400, height = 1600, units = 'px', dpi = 180)

# # epiR tool https://search.r-project.org/CRAN/refmans/epiR/html/epi.occc.html
# Div %>% select(KB, SM) %>% 
#   # rationale for scaling ?
#   mutate(across(everything(), ~scale(.x, scale = FALSE))) %>% 
#   epiR::epi.occc(pairs = TRUE)
# 
# # CCCRM tool
# Div_scaled_long <- Div %>% 
#   # Center only, because variances are equal
#   mutate(across(everything(), scale)) %>% 
#   rownames_to_column("Sample") %>% 
#   pivot_longer(values_to = divIDX, 
#                names_to = "tool",
#                cols = all_of(tools)) %>% 
#   mutate(tool = factor(tool, levels = tools)) %>% 
#   left_join(psSalivaKB@sam_data %>% data.frame %>% 
#               dplyr::select(treatDay) %>% 
#               rownames_to_column("Sample"),
#             by = "Sample")
# 
# ccc_result <- cccvc(Div_scaled_long %>% filter(tool !="MPA"), ry = divIDX, 
#                     rind = "Sample", rmet = "tool")
# summary(ccc_result)
# 
# ccc_result <- cccvc(Div_scaled_long %>% 
#                       filter(tool != 'KB'), 
#                     ry = "Shannon", 
#                     rind = "Sample", rmet = "tool")
# summary(ccc_result)
# 
# 
# 
# 
# 
# 
# ########################## Not updated below ###################################
# 
# # Dynamic formula : 
# formula <- reformulate('tool', response = divIDX)
# leveneTest(formula, data = Div_long) # same variance
# kruskal.test(formula, data = Div_long) # different distributions
# 
# #### Coefficient of variation ; testing whether they are different across categories
# # Coefficient of variation for each sample: 
# apply(Div, 2, function(x) (sd(x)/mean(x))) # Are they really different ?
# # Using Feltz and Miller 1996 implemented in : 
# # https://cran.r-project.org/web/packages/cvequality/vignettes/how_to_test_CVs.html
# mslr_test(1e5, Div_long[[divIDX]], Div_long$tool) # Same variation coefficient
# 
# #############################
# ### Bland & Altman Analysis ###
# #############################
# 
# BA_analysis <- compute_meandiff(Div, 'KB', 'SM')
# (test.norm <- shapiro.test(BA_analysis$Diff) %$% p.value)
# 
# # If differences are not normal, log transform diversity counts :
# if(test.norm <0.05) {
#   BA_analysis <- Div %>% 
#     mutate(across(everything(), log)) %>% 
#     compute_meandiff('KB', 'SM')
#   shapiro.test(BA_analysis$Diff) %$% p.value
# }
# 
# # global statistics of means and differences:
# mean_diff <- mean(BA_analysis$Diff)
# s <- sd(BA_analysis$Diff)
# upLimit = mean_diff + 1.96*s
# loLimit = mean_diff - 1.96*s
# n = length(BA_analysis$Diff)
# mean_sd = sqrt(s^2/n)
# s_sd = sqrt((3*s^2)/n)
# 
# # proportion of differences outside limits
# mean(BA_analysis$Diff < loLimit | BA_analysis$Diff > upLimit) 
# 
# # Bland-Altman Plot :
# ggplot(BA_analysis) +
#   theme_light() +
#   annotate(geom = 'rect', xmin = -Inf, xmax = Inf, ymin = upLimit - s_sd, 
#            ymax =  upLimit + s_sd, alpha = 0.2, fill = 'blue') +
#   annotate(geom = 'rect', xmin = -Inf, xmax = Inf, ymin = mean_diff - mean_sd, 
#            ymax = mean_diff + mean_sd, alpha = 0.2, fill = 'red') +
#   annotate(geom = 'rect', xmin = -Inf, xmax = Inf, ymin = loLimit - s_sd, 
#            ymax = loLimit + s_sd, alpha = 0.2, fill = 'blue') +
#   geom_hline(yintercept = mean_diff) +
#   geom_hline(yintercept = upLimit, linetype = 'dashed') +
#   geom_hline(yintercept = loLimit, linetype = 'dashed') +
#   geom_point(aes(x = mean, y = Diff)) +
#   labs(x = "Mean", y = "Difference", title = "Limits of Agreement between Kraken and Sourmash")
# 
# #### Trying out something... plot 3-way variation against mean
# sd_mean <- Div %>% rowwise %>% 
#   transmute(div_mean = mean(c_across(everything())),
#             div_sd = sd(c_across(everything()))) 
# 
# cor.test(sd_mean$div_mean, sd_mean$div_sd, tool = 'spearman') # no correlation
# plot(div_sd~div_mean, data = sd_mean)
# 
# ####################################################
# ### (Overall) Concordance Correlation Coefficient ###
# # Should we scale the diversities for the CCC test ? Only center?
# ####################################################
# 
# # epiR tool https://search.r-project.org/CRAN/refmans/epiR/html/epi.occc.html
# Div %>% select(KB, SM) %>% 
#   # rationale for scaling ?
#   mutate(across(everything(), ~scale(.x, scale = FALSE))) %>% 
#   epiR::epi.occc(pairs = TRUE)
# 
# # CCCRM tool
# Div_scaled_long <- Div %>% 
#   # Center only, because variances are equal
#   mutate(across(everything(), scale)) %>% 
#   rownames_to_column("Sample") %>% 
#   pivot_longer(values_to = divIDX, 
#                names_to = "tool",
#                cols = all_of(tools)) %>% 
#   mutate(tool = factor(tool, levels = tools)) %>% 
#   left_join(psSalivaKB@sam_data %>% data.frame %>% 
#               dplyr::select(treatDay) %>% 
#               rownames_to_column("Sample"),
#             by = "Sample")
# 
# ccc_result <- cccvc(Div_scaled_long %>% filter(tool !="MPA"), ry = divIDX, 
#                     rind = "Sample", rmet = "tool")
# summary(ccc_result)
# 
# ccc_result <- cccvc(Div_scaled_long %>% 
#                       filter(tool != 'KB'), 
#                     ry = "Shannon", 
#                     rind = "Sample", rmet = "tool")
# summary(ccc_result)
