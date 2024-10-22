library(pacman)
p_load(magrittr, tidyverse, 
       RColorBrewer, ggbeeswarm2, UpSetR, patchwork)

source('scripts/myFunctions.R')
source('scripts/5_DAA_fun.R')

# Full data
DAA_files <- 'Out/DAA*/*.tsv'
DAA <- Sys.glob(DAA_files) %>% map(read_tsv) %>% list_rbind

# NAFLD subset of significant taxa
subset_DAA_05_NAFLD_F <- DAA %>% 
   filter(taxRank == 'Family' &
      dataset == 'NAFLD' &
      adj.p <0.05) 

######################
### Taxon PA plot by DAA/db combination
######################

# tool_names <- subset_toolpairs %>% arrange(desc(DAA_tool)) %$% DAA_tool %>% unique
tool_names <- c(#"radEmu", 
                'corncob','MaAsLin2', 'Aldex2', 'ANCOMBC2', 'DESeq2', 'edgeR')
#db_names <- subset_toolpairs %>%  arrange(desc(database)) %$% database %>% unique
db_names <- c("KB20", "SM_gtdb-rs214-rep", "SM_genbank-2022.03", 
              "MPA_db2023", "MPA_db2019" ,"MOTUS" )    

# Filter out some tools to alleviate the plot
subset_toolpairs <- subset_DAA_05_NAFLD_F %>% 
  filter(DAA_tool %in% tool_names &
           database %in% db_names)

# Top taxa presence across combinations
top_taxa <- subset_toolpairs %>% 
  group_by(Taxon) %>% 
  dplyr::summarise(n = n()) %>% 
  filter(!str_detect(Taxon, '\\[')) %>% 
  arrange(desc(n)) %>%
  head(n = 20) %>% pull(Taxon)

# Data for the bottom matrix, two lines per pair string
tool_pairs <- expand_grid(tool_names, db_names) %>% 
  mutate(pair = paste(tool_names, db_names, sep= '_')) %>% 
  pivot_longer(cols = c('db_names', 'tool_names'), values_to = 'y') %>% 
  mutate(y = factor(y, levels = c(db_names, tool_names))) 

# Data for top plot matrix, 
tax_tool_pairs <- subset_toolpairs %>% 
  filter(Taxon %in% top_taxa) %>% 
  mutate(pair = paste(DAA_tool, database, sep= '_')) %>% 
  # Coefficients sign
  dplyr::select(Taxon, pair, coef) %>% 
  mutate(coef = factor(case_when(coef>0 ~ 1, coef<0 ~ -1)))

# Factor levels : taxa sorted by overall count
tax_count <- tax_tool_pairs %>% 
  dplyr::count(Taxon) %>% 
  arrange(desc(n)) %>% 
  pull(Taxon) 

# List all combinations
# comb_count <- tax_tool_pairs %>%
#   dplyr::count(pair, sort = TRUE) %>% 
#   pull(pair)

comb_factor <- tool_pairs %$% pair %>% unique %>% rev

# Refactor top matrix
tax_tool_pairs %<>% 
  mutate(pair = factor(pair, levels = comb_factor),
         Taxon = factor(Taxon, levels = tax_count))

# Refactor tool pairs
CCE_names_noreturn <- gsub("\n", " ", CCE_names)
tool_pairs %<>%
  mutate(pair = factor(pair, levels = comb_factor)) %>% 
  filter(pair %in% tax_tool_pairs$pair) %>% 
  # Rename labels
  mutate(y = recode(y, !!!CCE_names_noreturn))

# Plot !
top_plot <- tax_tool_pairs %>% 
  ggplot(aes(x = pair, y = Taxon)) +
  geom_point(aes(colour = coef), size = 3) +
  scale_colour_manual(values = c('red3', "blue3"),
                      labels = c('Negative', 'Positive')) +
  theme_minimal() +
  labs(colour = 'Taxon association')

bottom_plot <- tool_pairs %>% 
  ggplot(aes(x = pair, y = y)) +
  geom_hline(yintercept = match("radEmu", levels(tool_pairs$y)) - 0.5, 
             color = "black", linewidth = 0.2, linetype = 'solid') +
  geom_line(aes(group = pair)) +
  geom_point(aes(colour = name), size = 3) +
  scale_colour_manual(values = c('orange2', 'purple2'), 
                      labels = c('Community composition estimation', 
                                 'Differential abundance'))+
  labs(colour = 'Tool type') +
  theme_minimal() 

(top_plot / bottom_plot) +
  plot_layout(guides = "collect",
              heights = c(length(tax_count),
                          length(c(tool_names, db_names)))) &
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title.position = 'top') 

ggsave('Out/CSHL_poster/DAA_Story.pdf', bg = 'white', width = 2000, height = 1600, units = 'px', dpi = 180)

#################
## Clustering ####
###################

# Add unique tool-set variable
subset_DAA_05_NAFLD_F %<>% 
  dplyr::mutate(present = 1,
                tool_set = paste0(DAA_tool, '__', database)) 


# Some contradictions exist (taxa as DAA with different signs) so 
# we remove those taxa when they exist (will underestimate dissimilarity)
subset_no_contradiction <- subset_DAA_05_NAFLD_F %>% 
  group_by(Taxon) %>% 
  filter(!(any(coef > 0) & any(coef < 0))) %>% 
  ungroup()
  
  
# PA matrix
PA <- subset_no_contradiction %>% 
  filter(DAA_tool %in% tool_names ) %>% 
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
   # stat_ellipse(level=0.85, geom = "polygon", alpha = 0.1, 
   #              aes(colour = CCE_approach), fill = NA) + 
  ggbeeswarm2::geom_beeswarm(aes(colour = database, shape = DAA_tool), 
                size = 3, stroke = 1.5, method = 'swarm2', spacing = 1.5) +
  scale_shape_manual(values = setNames(DAA_metadata$plot_shape, DAA_metadata$DAA_tool)) +
  scale_colour_manual(values = tool_colours, labels = CCE_names) +
  theme_minimal() +
  # geom_text(data = approach_text, aes(x = X, y = Y, label = CCE_approach), 
  #           color = "black", size = 5, fontface = "bold") +
  labs(x = paste0('PCo1 (', eig[1],')'), 
       y = paste0('PCo1 (', eig[2],')'))

ggsave('Out/CSHL_poster/pcoa_tool_sets.pdf', bg = 'white', width = 2400, height = 2000, units = 'px', dpi = 180)

# Permanova?
adonis2(formula = dist ~ taxonomy + CCE_approach + Taxon_bias + Compositional, 
        permutations = 999,
        data = tool_sets_meta,
        by = 'margin',
      #  na.action = na.exclude,
        parallel = 8)

  
# Long DF with every possible tool combination for each taxon
# 
# # P/A matrix for datadases 
# wide_db <- subset_DAA_05_NAFLD_F %>% 
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
#   filter(sum(c_across(all_of(tool_names))) == 1 &
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


