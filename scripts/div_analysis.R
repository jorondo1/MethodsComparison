library(pacman)
p_load(phyloseq, tidyverse,cccrm, magrittr, cvequality)
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

Div_long <- bind_rows(Species = compile_diversity(ps_species.ls), 
          Genus = compile_diversity(ps_genus.ls), 
          Family = compile_diversity(ps_family.ls),
          .id = 'Rank')

# Visualise differences in diversity across tools
Div_long %>% 
  filter(Rank != 'Family' & index == 'H_0' & dataset != 'Feces') %>% 
  ggplot(aes(x = database, y = value, fill = database)) +
  geom_boxplot(outlier.size = 0.5, size = 0.3) + geom_line(aes(group = Sample), alpha=0.3, linewidth = 0.2)+ theme_light() +
  facet_grid(cols = vars(dataset), rows=vars(index, Rank), scales = 'free') +
  theme(axis.text.x = element_blank())

ggsave('Out/diversity_filt.pdf', bg = 'white', 
       width = 1900, height = 2400, units = 'px', dpi = 180)

# Check distribution of indices :
# lapply(div_rare, function(sublist) {
#   lapply(sublist, function(subsublist) {
#     lapply(subsublist, function(element) {
#       shapiro.test(element)$p.value
#     })})
# }) 

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

# 
# Function to apply cccvc and return results in a tibble
cccvc_groups <- function(df) { # run cccvc on filtered data: 
  ccc_out <- cccvc(df, ry = 'value', rind = 'Sample', rmet = 'database')  
  # compile
  tibble(CCC = ccc_out$ccc[1],      # CCC 
        LL_CI_95 = ccc_out$ccc[2], # lower limit CI
        UL_CI_95 = ccc_out$ccc[3], # upper limit CI
        SE_CCC = ccc_out$ccc[4])   # standard error
}

# Apply across all dataset & index combinations
ccc.df <- Div %>%
  dplyr::filter(database %in% c("MPA_db2023", "KB51", "MOTUS", "SM_gtdb_rs214_full")) %>% 
  group_by(dataset, index, Rank) %>%         # group data by dataset and index
  group_modify(~ cccvc_groups(.x)) %>% # apply cccvc_groups to each group
  ungroup                              # remove grouping structure

ggplot(ccc.df, aes(x = factor(index), y = CCC, color = Rank)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +  # Dodge points
  geom_errorbar(aes(ymin = LL_CI_95, ymax = UL_CI_95), 
                width = 0.2, position = position_dodge(width = 0.5)) +  # Dodge error bars
  facet_wrap(~ dataset) +  # Facet by dataset
  labs(x = "Index", y = "CCC", title = "Concordance Correlation Coefficient") +  # Axis and title labels
  theme_minimal() + ylim(c(0,1))

# Positive control
Div %>%
  dplyr::filter(database %in% c("SM_gtdb_rs214_full", "SM_genbank_202203") &
                  Rank == 'Species') %>% 
  group_by(dataset, index) %>%         # group data by dataset and index
  group_modify(~ cccvc_groups(.x)) %>% # apply cccvc_groups to each group
  ungroup %>%                          # remove grouping structure

ggplot(aes(x = factor(index), y = CCC, color = index)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LL_CI_95, ymax = UL_CI_95), width = 0.2) +  # Confidence interval
  facet_wrap(~ dataset) +  # Facet by dataset
  labs(x = "Index", y = "CCC", title = "Concordance Correlation Coefficient") +  # Axis and title labels
  theme_minimal() + ylim(c(0,1))


################################
### CCC across pairs of tools ##
##################################

Div_saliva <- Div %>% filter(dataset == 'Saliva')

# Function to apply cccvc to a pair of tools, with error handling
cccvc_pairwise <- function(df, tool_pair) {
  # Filter the data to keep only the two tools in the pair
  df_pair <- df %>% filter(database %in% tool_pair)
  
  # Try to compute cccvc, and return NA in case of an error
  tryCatch({
    # Run cccvc and return results in a tibble
    ccc_out <- cccvc(df_pair, ry = 'value', rind = 'Sample', rmet = 'database')
    
    tibble(
      tool1 = tool_pair[1],       # Correct tool names
      tool2 = tool_pair[2], 
      CCC = ccc_out$ccc[1],       # CCC 
      LL_CI_95 = ccc_out$ccc[2],  # Lower limit CI
      UL_CI_95 = ccc_out$ccc[3],  # Upper limit CI
      SE_CCC = ccc_out$ccc[4]     # Standard error
    )
  }, error = function(e) {
    # If there is an error, return NA values for this pair
    tibble(
      tool1 = tool_pair[1],
      tool2 = tool_pair[2],
      CCC = NA,           # NA for the CCC
      LL_CI_95 = NA,      # NA for lower limit CI
      UL_CI_95 = NA,      # NA for upper limit CI
      SE_CCC = NA         # NA for standard error
    )
  })
}

# Apply the pairwise CCC calculation for each dataset and index group
ccc_pairwise_df <- Div_saliva %>%
  group_by(dataset, index, Rank) %>%
  group_modify(~ {
    # Extract the current group
    current_group <- .x
    
    # Create unique tool pairs within the group
    tool_pairs <- unique(current_group$database) %>% combn(2, simplify = FALSE)
    
    # For each pair, apply the cccvc_pairwise function
    map_dfr(tool_pairs, function(pair) cccvc_pairwise(current_group, pair))
  }) %>%
  ungroup()

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
