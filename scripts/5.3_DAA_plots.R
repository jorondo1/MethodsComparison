library(pacman)
p_load(magrittr, tidyverse, 
       RColorBrewer, ggbeeswarm2, UpSetR, patchwork)

source('scripts/myFunctions.R')
source('scripts/5_DAA_fun.R')

# Full data
DAA <- rbind(
  read_tsv('Out/DAA/Maaslin2.tsv'),
  read_tsv('Out/DAA/AncomBC2.tsv'),
  read_tsv('Out/DAA/edgeR.tsv'), # Too many taxa, needs dealing with !
  read_tsv('Out/DAA/DESEq2.tsv'),
  read_tsv('Out/DAA/radEmu.tsv'),
  read_tsv('Out/DAA/Aldex2.tsv'),
  read_tsv('Out/DAA/ZicoSeq.tsv')
)

# NAFLD subset of significant taxa
subset_DAA_05 <- DAA %>% 
   filter(taxRank == 'Family' &
      dataset == 'NAFLD' &
      adj.p <0.05) 

### Wide matrix for upset plot 
wide_DAA <- subset_DAA_05 %>% 
  dplyr::select(Taxon, DAA_tool) %>% 
  distinct() %>% 
  mutate(present = 1) %>% 
  pivot_wider(names_from = DAA_tool, values_from = present, values_fill = 0) 

# Simple upset plot 
upset(wide_DAA %>% as.data.frame, 
      sets = names(wide_DAA)[-1],
      order.by = "freq")

######################
### Taxon PA plot by DAA/db combination
######################

# Filter out some tools to alleviate the plot
subset_toolpairs <- subset_DAA_05 %>% 
  filter(!DAA_tool %in% c('edgeR','ZicoSeq') &
           !database %in% c('MPA_db2022','KB51'))

# Top taxa presence across combinations
top_taxa <- subset_toolpairs %>% 
  group_by(Taxon) %>% 
  dplyr::summarise(n = n()) %>% 
  filter(!str_detect(Taxon, '\\[')) %>% 
  arrange(desc(n)) %>%
  head(n = 20) %>% pull(Taxon)

# Extract tool names
tool_names <- subset_toolpairs %>% arrange(desc(DAA_tool)) %$% DAA_tool %>% unique
db_names <- subset_toolpairs %>% arrange(desc(database)) %$% database %>% unique

# Data for the bottom matrix, two lines per pair string
tool_pairs <- expand_grid(tool_names, db_names) %>% 
  mutate(pair = paste(db_names, tool_names, sep= '_')) %>% 
  pivot_longer(cols = c('db_names', 'tool_names'), values_to = 'y') %>% 
  mutate(y = factor(y, levels = c(db_names, tool_names))) 

# Data for top plot matrix, 
tax_tool_pairs <- subset_toolpairs %>% 
  filter(Taxon %in% top_taxa) %>% 
  mutate(pair = paste(database, DAA_tool, sep= '_')) %>% 
  # Coefficients sign
  dplyr::select(Taxon, pair, coef) %>% 
  mutate(coef = factor(case_when(coef>0 ~ 1, coef<0 ~ -1)))

# List all combinations
comb_count <- tax_tool_pairs %>%
  dplyr::count(pair, sort = TRUE) %>% 
  pull(pair)

# Factor levels : taxa sorted by overall count
tax_count <- tax_tool_pairs %>% 
  dplyr::count(Taxon) %>% 
  arrange(desc(Taxon)) %>% 
  pull(Taxon) 

# Refactor top matrix
tax_tool_pairs %<>% 
  mutate(#pair = factor(pair, levels = comb_count),
         Taxon = factor(Taxon, levels = tax_count))

# Refactor tool pairs
tool_pairs %<>%
#  mutate(pair = factor(pair, levels = comb_count)) %>% 
  filter(pair %in% tax_tool_pairs$pair) %>% 
  # Rename labels
  mutate(y = recode(y, !!!CCE_names))

# Plot !
top_plot <- tax_tool_pairs %>% 
  ggplot(aes(x = pair, y = Taxon)) +
  geom_point(aes(colour = coef), size = 3) +
  scale_colour_manual(values = c('#9E0142', "#5E4FA2"),
                      labels = c('Negative', 'Positive')) +
  theme_minimal() +
  labs(colour = 'Taxon association')

bottom_plot <- tool_pairs %>% 
  ggplot(aes(x = pair, y = y)) +
  geom_hline(yintercept = match("radEmu", levels(tool_pairs$y)) - 0.5, 
             color = "black", linewidth = 0.2, linetype = 'solid') +
  geom_line(aes(group = pair)) +
  geom_point(aes(colour = name), size = 3) +
  scale_colour_manual(values = c('#FDAE61', '#66C2A5'), 
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

ggsave('Out/Story.pdf', bg = 'white', width = 2000, height = 1600, units = 'px', dpi = 180)

# Long DF with every possible tool combination for each taxon
# 
# # P/A matrix for datadases 
# wide_db <- subset_DAA_05 %>% 
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


