library(pacman)
p_load(magrittr, tidyverse, 
       RColorBrewer, ggbeeswarm2, UpSetR, patchwork, cowplot)

source('scripts/myFunctions.R')

###########################################################
#
#  /$$$$$$$  /$$        /$$$$$$  /$$$$$$$$         /$$  
# | $$__  $$| $$       /$$__  $$|__  $$__/       /$$$$  
# | $$  \ $$| $$      | $$  \ $$   | $$         |_  $$  
# | $$$$$$$/| $$      | $$  | $$   | $$           | $$  
# | $$____/ | $$      | $$  | $$   | $$           | $$  
# | $$      | $$      | $$  | $$   | $$           | $$  
# | $$      | $$$$$$$$|  $$$$$$/   | $$          /$$$$$$
# |__/      |________/ \______/    |__/         |______/
#
# Taxon upset plot by DAA/db combination
#
###########################################################

######################
### 1. Subset data ####
########################

if(!file.exists('Out/_Rdata/taxa_relAb_metrics_flat.RDS')) {
  source('scripts/5.0_taxa_relAb_metrics_flat.R')
}
taxa_relAb_metrics <- readRDS('Out/_Rdata/taxa_relAb_metrics.RDS')

# Import data
DAA <- Sys.glob('Out/DAA/*.tsv') %>% 
  map(read_tsv, show_col_types = FALSE) %>% 
  list_rbind()

# Data to subset
which_tools <- c('DESeq2', 'edgeR', 'ZicoSeq', 'ANCOMBC2', 'radEMU', 'Aldex2','corncob','MaAsLin2')
which_databases <- c('KB90','KB45',  'SM_RefSeq_20250528', 
                     'SM_gtdb-rs220-rep', 'KB45_GTDB', 'KB90_GTDB', 
                      'MPA_db2023', 'MOTUS')

# Subset to dataset + Taxon
which_taxrank <- 'Genus'
which_dataset <- 'PD'
taxa_relAb_metrics_flat <- taxa_relAb_metrics[[which_taxrank]][[which_dataset]]

DAA_subset <- DAA %>% 
  filter(taxRank == which_taxrank &
           dataset == which_dataset) %>% 
  select(-taxRank, -dataset)

# Filter out some tools to alleviate the plot
subset_toolpairs <- DAA_subset %>% 
  filter(DAA_tool %in% which_tools 
         & database %in% which_databases
  )

######################
### 2. Select taxa ####
######################## 
# We want the most overall abundant taxa (across all CCE method) that have
# been called significant by at least one tool. 

# Taxa found significant in at least N methodologies
p_threshold <- 0.01
taxa_sig <- subset_toolpairs %>% 
  filter(adj.p <= p_threshold) %>% 
  group_by(Taxon) %>% 
  summarise(n = n()) %>% 
  filter(n >= 3) %>% 
  pull(Taxon) %>% unique()

# List taxa found by every CCE method
taxa_in_all_db <- taxa_relAb_metrics_flat %>% 
  filter(database %in% which_databases
         & Taxon %in% taxa_sig) %>% 
  group_by(Taxon) %>% 
  dplyr::summarise(n = n()) %>% 
  arrange(desc(n)) %>% 
  filter(n == length(which_databases)) %>% 
#  filter(n >=2) %>% 
  pull(Taxon)

# List taxa by top overall mean relative abundance
top_taxa <-  taxa_relAb_metrics_flat %>% 
  filter(Taxon %in% taxa_in_all_db) %>% # found by all CCE methods
  group_by(Taxon) %>% 
  summarise(overall_meanRelAb = mean(meanRelAb)) %>% 
  arrange(desc(overall_meanRelAb)) %>% 
 # head(n = 30) %>% 
  pull(Taxon)

#################################
### 3. Create plot dataframes ####
###################################
# For each database*DAA_tool combination, create a unique identifier
# with both variables in other columns; the pair will serve as the x axis.
# Then, generate the data for the main plot, with coefficients converted
# to negative/positive binary.

# Data for the bottom matrix, two lines per pair string
tool_pairs <- expand_grid(x = which_tools, y = which_databases) %>% 
  mutate(pair = paste(x, y, sep= '_')) %>% 
  mutate(y = factor(y, levels = which_databases),
         x = factor(x, levels = which_tools))

# Tool pair factors in order
comb_factor <- rev(tool_pairs$pair) # %>% droplevels

# Data for main plot matrix
tax_tool_pairs <- subset_toolpairs %>% 
  filter(Taxon %in% top_taxa) %>% 
  mutate(pair = paste(DAA_tool, database, sep= '_')) %>% 
  # Convert coefficients to signs: 
  mutate(coef = factor(case_when(coef>0 ~ 1, coef<0 ~ -1)))

# Factor levels : taxa sorted by overall agreement score
tax_levels <- tax_tool_pairs %>% 
  group_by(Taxon) %>% 
  dplyr::summarise(score = abs(sum(as.numeric(as.character(coef))))) %>%
  arrange(desc(score)) %>% 
  pull(Taxon)

# Refactor tool pairs
CCE_names_noreturn <- gsub("\n", " ", CCE_names) # remove returns from full tool names 
tool_pairs %<>%
  mutate(pair = factor(pair, levels = comb_factor)) %>% 
  dplyr::filter(pair %in% tax_tool_pairs$pair) %>% 
  # Rename labels
  mutate(y = recode(y, !!!CCE_names_noreturn))

# Build the real plot data tibble
main_plot_data <- taxa_relAb_metrics_flat %>%  # Add the scaled_diff and meanRelAb variables
  dplyr::filter(Taxon %in% tax_levels) %>% 
  left_join(tax_tool_pairs, ., 
            by = c('Taxon', 'database')) %>% 
  dplyr::filter(adj.p < p_threshold) %>% 
  mutate(pair = factor(pair, levels = comb_factor),
         Taxon = factor(Taxon, levels = tax_levels),
         DAA_tool = factor(DAA_tool, levels = which_tools))

###############
### 4. PLOT ####
#################

# Dataframe-generating function for a striped background with facets
set_background_alpha <- 0.4
striped_background_facet <- function(data_pairs, x_var, facet_var) {
  background_data <- data_pairs %>%
    distinct(.data[[facet_var]], .data[[x_var]], database) %>%
    left_join(CCE_metadata %>% select(Database, plot_colour),
              by = c('database' = 'Database')) %>%
    group_by(.data[[facet_var]]) %>%
    mutate(
      x_index = match(database, rev(which_databases))
    ) %>%
    ungroup() 

  #Return a list of ggplot layers
  list(
    geom_rect(
      data = background_data,
      aes(
        xmin = x_index - 0.5,
        xmax = x_index + 0.5,
        ymin = -Inf,
        ymax = Inf,
        fill = plot_colour
      ),
      inherit.aes = FALSE,
      alpha = set_background_alpha
    ),
    scale_fill_identity()
  )
}

# PLOT !

main_plot <- main_plot_data %>% 
  ggplot(aes(x = pair, y = Taxon)) +
  striped_background_facet(main_plot_data, 'pair', 'DAA_tool') +
  geom_point(aes(shape = coef), size = 1.8, colour = 'grey10') +
  facet_grid(.~DAA_tool, scales = 'free', space = 'free')+
  scale_shape_manual(values = c(1,19),
                      labels = c('Group A', 'Group B')) +
  #scale_fill_manual(values = c('grey','black'))+
  #scale_size_continuous(breaks = c(0.2, 0.7), range = c(0.2,4)) +
  theme_light() +
  labs(shape = 'Taxon association', 
       size = 'Scaled difference in mean abundance between groups:') +
  theme(plot.margin = margin(t = 0, unit = "cm") ,
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = 'grey50'),
        plot.title = element_text(margin = margin(b = 0)) ,
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.justification = 'left',
        legend.title.position = 'top',
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.spacing.x = unit(0, "lines"),
        legend.margin = margin(0,0,0,0)
       # panel.border = element_blank()
  ); main_plot

ggsave(plot = main_plot, 'Out/ISMB2025/DAA_Story_Genus.pdf', 
       bg = 'white', width = 2000, height = 1400, units = 'px', dpi = 180)

# Dummy legend
plot_colours_df <- CCE_metadata %>% 
  filter(Database %in% which_databases) %>% 
  select(Database, plot_colour, MethodName) %>% 
  mutate(Database = factor(Database, levels = rev(which_databases)))

db_colours <- plot_colours_df %>% select(Database, plot_colour) %>% deframe()

p_legend_dummy <- plot_colours_df %>% 
  ggplot(aes(x = 1, y = 1, fill = Database)) +
  geom_tile(width = 0.1, height = 0.1, alpha = set_background_alpha) + # A tiny geom to trick ggplot into making a legend
  scale_fill_manual(
    name = 'Community Composition Estimation Methodology',
    values = db_colours,
    labels = CCE_names
  ) +
  theme_void() + # Remove all plot elements (axes, titles, background, etc.)
  theme(
    legend.position = 'bottom', # Position for the extracted legend
    legend.direction = "horizontal",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15, face = "bold"),
    legend.title.position = 'top',
    legend.key.size = unit(0.7, "cm") # Adjust the size of the legend keys
  ) +
  guides(fill = guide_legend(nrow = 2)) # Ensure legend squares are opaque, not alpha=0.5

ggsave(plot = p_legend_dummy, 'Out/ISMB2025/DAA_Story_Genus_legend.pdf', 
       bg = 'white', width = 1550, height = 170, units = 'px', dpi = 180)










#################
## Clustering ####
###################

# Add unique tool-set variable
DAA_subset %<>% 
  dplyr::mutate(present = 1,
                tool_set = paste0(DAA_tool, '__', database)) 

# Some contradictions exist (taxa as DAA with different signs) so 
# we remove those taxa when they exist (will underestimate dissimilarity)
subset_no_contradiction <- DAA_subset %>% 
  group_by(Taxon) %>% 
  filter(!(any(coef > 0) & any(coef < 0))) %>% 
  ungroup()

# PA matrix
PA <- subset_no_contradiction %>% 
  filter(DAA_tool %in% which_tools ) %>% 
  dplyr::select(Taxon, tool_set, present) %>% 
  pivot_wider(names_from = tool_set, values_from = present, values_fill = 0) %>% 
  column_to_rownames('Taxon') %>% t 
###### Subset to taxa common across db ??

# Bray-curtis on PA = Sørenson dice (more weight to common)
dist <- PA %>% vegdist(method = 'bray')
pcoa_res <- capscale(dist~1, distance = 'bray')

eig <- (round(pcoa_res$CA$eig[1:3]/sum(pcoa_res$CA$eig),3)*100) %>% paste0('%')
message(paste("First 3 PCo :",eig[1], ',', eig[2], ',', eig[3]))

# create output list
tool_sets_meta <- subset_no_contradiction %>% 
  dplyr::select(tool_set, DAA_tool, database) %>% 
  distinct %>% 
  inner_join(scores(pcoa_res)$sites[,1:2] %>% 
               as.data.frame %>% 
               rownames_to_column('tool_set'),
             by = 'tool_set') %>% 
  left_join(CCE_metadata, by = 'database') %>% 
  left_join(DAA_metadata, by = 'DAA_tool') %>% 
  mutate(DAA_tool = factor(DAA_tool, levels = names(tool_vars)),
         database = factor(database, levels = names(tool_colours)))

approach_text <- tibble(
  CCE_approach = c('DNA-to-Marker', 'DNA-to-DNA'),
  X = c(0, -1),
  Y = c(1.2, -1.2)
)

tool_sets_meta %>% 
  ggplot(aes(x = MDS1, y = MDS2)) +
  # stat_ellipse(level=0.7, geom = "polygon", alpha = 0.1, aes(colour = CCE_approach), fill = NA) + 
  ggbeeswarm2::geom_beeswarm(aes(colour = database, shape = DAA_tool), 
                             size = 4, stroke = 1.5, method = 'swarm2', spacing = 1.5) +
  scale_shape_manual(values = setNames(DAA_metadata$plot_shape, DAA_metadata$DAA_tool)) +
  scale_colour_manual(values = tool_colours, labels = CCE_names) +
  theme_light() +
  # geom_text(data = approach_text, aes(x = X, y = Y, label = CCE_approach), color = "black", size = 5, fontface = "bold") +
  labs(x = paste0('PCo1 (', eig[1],')'), 
       y = paste0('PCo2 (', eig[2],')'),
       shape = 'DAA tool', colour = 'Composition\nestimation tool') +
  guides(
    shape = guide_legend(order = 2),  # Put shape legend first
    colour = guide_legend(order = 1)    # Put fill legend second
  ) +
  theme(
    legend.position = c(0.016, 0.72),
    legend.justification = "left",                # Align the legend to the left
    legend.background = element_rect(fill = "white", color = "black", size = 0.5)
  )


ggsave('Out/CSHL_poster/pcoa_tool_sets.pdf', bg = 'white',
       width = 2400, height = 2400, units = 'px', dpi = 260)

# Permanova?
adonis2(formula = dist ~ taxonomy + CCE_approach + Taxon_bias, 
        permutations = 999,
        data = tool_sets_meta,
        by = 'margin',
        #  na.action = na.exclude,
        parallel = 8)


# Long DF with every possible tool combination for each taxon
# 
# # P/A matrix for datadases 
# wide_db <- DAA_subset %>% 
#   dplyr::select(Taxon, database) %>% 
#   distinct() %>% 
#   mutate(present = 1) %>% 
#   pivot_wider(names_from = database, values_from = present, values_fill = 0) 
# 
# # Merge w/ DAAtool mx
# wide_full <- full_join(wide_db, wide_DAA, by = 'Taxon')
# wide_full$comb <- apply(wide_full[,-1], 1, paste, collapse = "")
# 
# wide_pairs <- rowwise(wide_full) %>%
#   filter(sum(c_across(all_of(which_tools))) == 1 &
#            sum(c_across(all_of(db_names))) == 1) %>%
#   ungroup()
# 
# 
# wide_full %>% 
# find_strings_with_two_ones('comb', 9) %>% View
# 
# 
# 
# #######################
# ### UpsetUpset Plot ####
# #########################
# # Generate multiple tool combinations 
# create_combinations <- function(df) {
#   # Get numeric column names
#   col_names <- df %>% select(-Taxon) %>% colnames
#   # Generate all combinations of 2 or more columns
#   combinations <- map(2:length(col_names), 
#                       ~combn(col_names, .x, simplify = FALSE)) %>%
#     flatten()
#   # For each combination, create a new column with the product and a corresponding name
#   for (combo in combinations) {
#     combo_name <- paste(combo, collapse = "_")
#     df[[combo_name]] <- purrr::reduce(select(df, all_of(combo)), `*`)
#   }
#   return(df)
# }
# 
# # Long DF with every possible tool combination for each taxon
# DAA_comb_long <- wide_DAA %>% 
#   as_tibble %>% 
#   filter(Taxon %in% top_taxa) %>% 
#   create_combinations() %>% # Custom function
#   pivot_longer(cols = -Taxon, 
#                names_to = "toolComb", 
#                values_to = "Value") %>% 
#   dplyr::filter(Value == 1) %>% dplyr::select(-Value) 
# 
# 
# # DAA_comb_long %>%
# #   tidyr::separate(toolComb, into = paste0("Tool", 1:7), sep = "_", fill = "right", remove = FALSE) %>%
# #   pivot_longer(cols = starts_with("Tool"), names_to = "Tool_order", values_to = "Tool") %>%
# #   filter(!is.na(Tool)) %>% 
# #   # Plot ! 
# #   ggplot(aes(x = Tool_order, y = Tool, fill = Tool)) +
# #   geom_point(shape = 21, size = 5, color = "black") +
# #   scale_fill_manual(values = c("white", "black")) +
# #   labs(x = "Tool Combinations", y = "Tools") +
# #   theme_minimal() +
# #   theme(axis.text.x = element_text(angle = 90, hjust = 1),
# #         axis.title.x = element_blank())  # Keep the x-axis labels consistent with the top plot
# # 
#   
# # BOTTOM starting from upsetPlot:
# upset_data <- UpSetR::upset(
#   wide_DAA, 
#   sets = names(wide_DAA)[-1],  # Exclude Taxon column
#   order.by = "freq",
#   nsets = length(names(wide_DAA)[-1])# Number of set
# )
# 
# upset_comb <- upset_data$New_data %>%
#   mutate(toolComb = apply(upset_data$New_data[ , -1], 1, function(x) paste(names(x)[which(x == 1)], collapse = "_"))) %>%
#   distinct(Taxon, toolComb) %>%
#   dplyr::filter(Taxon %in% top_taxa)
# 
# upset_comb%>% 
# ggplot(aes(x = toolComb, y = Taxon)) +
#   geom_point() +
#   labs(x = "Tool Combinations", y = "Taxa") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))


