library(pacman)
p_load(magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
       tidyverse, kableExtra,
       rstatix)

ps_rare.ls <- read_rds('Out/_Rdata/ps_rare.ls.rds')
source("scripts/0_Config.R")
theme_set(theme_light())


#   $$\ $$\         $$$$$$$$\ $$$$$$\  $$$$$$\             $$$$$$\  
#   $$ \$$ \        $$  _____|\_$$  _|$$  __$$\           $$  __$$\ 
# $$$$$$$$$$\       $$ |        $$ |  $$ /  \__|          \__/  $$ |
# \_$$  $$   |      $$$$$\      $$ |  $$ |$$$$\            $$$$$$  |
# $$$$$$$$$$\       $$  __|     $$ |  $$ |\_$$ |          $$  ____/ 
# \_$$  $$  _|      $$ |        $$ |  $$ |  $$ |          $$ |      
#   $$ |$$ |        $$ |      $$$$$$\ \$$$$$$  |$$\       $$$$$$$$\ 
#   \__|\__|        \__|      \______| \______/ \__|      \________|
#                         

###################
# Beta diversity ###
#####################


# Procruste analysis ------------------------------------------------------



# Another way of looking at it :  -----------------------------------------
# This is the first way I came up with, before being told that
# procrustes would be better



# using collapsed pairwise matrices (BC and rAitchison):
# Sample_pair, Dataset, idx, idx_value, CCE, Approach

these_datasets <- c('Moss', 'NAFLD', 'P19_Gut', 'P19_Saliva', 'PD', 'Bee', 'AD_Skin')

# Define list of pairs of interest, with names used for facet_grid
db_pairs_eval <- list(
  `A. Kraken 0.45 :\nGTDB 220 – RefSeq` = c('KB45_GTDB', 'KB45'),
  `B. Sourmash :\nGTDB 220 – RefSeq` = c('SM_gtdb-rs220-rep', 'SM_RefSeq_20250528'),
  `C. Sourmash GTDB 220 –\n Kraken RefSeq` = c('SM_gtdb-rs220-rep', 'KB45'),
  `D. Kraken GTDB 220 –\n Sourmash RefSeq` = c('KB45_GTDB', 'SM_RefSeq_20250528'),
  `E. GTDB 220 :\nSourmash – Kraken 0.45` = c('SM_gtdb-rs220-rep', 'KB45_GTDB'),
  `F. RefSeq :\nSourmash – Kraken 0.45` = c('SM_RefSeq_20250528', 'KB45'),
  `G. DNA-to-Marker :\nmOTUs3 – MetaPhlAn4` = c('MOTUS', 'MPA_db2023')
)

db_pairs_ctrl <- list(
  `A. Sourmash\nGTDB220 – GTDB214` = c('SM_gtdb-rs220-rep', 'SM_gtdb-rs214-rep'),
  `B. MetaPhlAn versions\n2023 – 2022` = c('MPA_db2023', 'MPA_db2022'),
  `C. GTDB taxonomy\n214 Full – 214 Reps (using Sourmash)` = c('SM_gtdb-rs214-full','SM_gtdb-rs214-rep'),
  `D. NCBI taxonomy\nGenbank – RefSeq (using Sourmash)` = c('SM_genbank-2022.03', 'SM_RefSeq_20250528')
)

# Load data 
pairwise_distances <- read_rds('Out/_Rdata/pairwise_dist.RDS') %>% 
  left_join(CCE_metadata, 
            by = 'Database') %>% 
  filter(Dist == 'bray' 
         #& Database %in% these_databases
         & Dataset %in% these_datasets
  ) %>% select(-Dist)

# Function to compute differences in pairwise distances between two tools 
compute_distance_differences <- function(df, tool1, tool2) {
  df %>%
    filter(Database %in% c(tool1, tool2)) %>%
    #Create sample pair ID with consistence and uniqueness
    mutate(Pair = ifelse(Sample1 < Sample2, 
                         paste(Sample1, Sample2, sep = "_"),
                         paste(Sample2, Sample1, sep = "_"))) %>%
    select(Database, Dataset, Pair, Distance) %>%
    # Pivot wide to manually compute difference
    pivot_wider(names_from = Database, values_from = Distance) %>%
    filter(complete.cases(.)) %>%
    mutate(abs_diff = abs(.[[tool1]] - .[[tool2]]),
           dist_diff = .[[tool1]] - .[[tool2]]) %>% # instead of !!sym() , thanks deepseek
    select(-all_of(c(tool1, tool2)))
}

# Iterate over pairs of interest
pairwise_dist_gap.df <- imap(
  c(db_pairs_eval,db_pairs_ctrl),
  function(tool_pair, pair_name){
    compute_distance_differences(pairwise_distances, 
                                 tool_pair[1], tool_pair[2]) %>% 
      mutate(
        Pair_name = pair_name
      )
  }) %>% list_rbind() %>% 
  mutate(
    Pair_name = factor(
      Pair_name, 
      levels = names(c(db_pairs_eval,db_pairs_ctrl))
      )
    )

# Filter pairs, set factor levels
pw_dist_gap_eval.df <- pairwise_dist_gap.df %>% 
  filter(Pair_name %in% names(db_pairs_eval)) %>% 
  mutate(Dataset = factor(Dataset, 
                          levels = c('P19_Saliva', 'P19_Gut', 'NAFLD',  'PD',  'AD_Skin','Moss','Bee')
  ))
# Same for controls: 
pw_dist_gap_ctrl.df <- pairwise_dist_gap.df %>% 
  filter(Pair_name %in% names(db_pairs_ctrl))

# Summary
pairwise_dist_summary <- 
  rbind(pw_dist_gap_eval.df, pw_dist_gap_ctrl.df) %>% 
  group_by(Pair_name, Dataset) %>% 
  summarise(#mean_diff = mean(dist_diff),
    #sd_diff = sd(dist_diff),
    median_diff = median(dist_diff),
    mad_diff = mad(dist_diff),
    # min_diff = min(dist_diff),
    max_diff = max(dist_diff),
    #  rcv_diff = mad_diff/median_diff,
    .groups = 'drop'
  ) %>% 
  mutate(across(where(is.numeric), ~round(.x, 3)))

# Output summary table
pairwise_dist_summary %>% 
  mutate(Pair_name = str_replace(Pair_name, "\n", " ")) %>% # prevents markdown pipes from being added
  kable(align = "l") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "responsive")) %>% 
  save_kable(paste0('Out/memoire/tables/beta_diffs.html'))

# PLOT Boxplot with differences between pairs of interest
plot_dist_gap <- function(df){
  
  ggplot(df, aes(x = Dataset, y = dist_diff, fill = Dataset)) +
    geom_hline(aes(yintercept = 0), 
               color = "black", linewidth = 0.5, linetype = 'dashed') +
    geom_violin(linewidth = 0.2, draw_quantiles = c(0.5)) +
    facet_grid(.~Pair_name, scale = 'free') +
    theme_light() +
    scale_fill_manual(values= c("#b86092",  "#a40000", "#16317d", "#de722a", "#00b7a7", "#007e2f", "#ffcd12"),
                      labels = dataset_names)+ 
    labs(y = 'Same-pair differences in dissimilarities',
         fill = 'Dataset') +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = 'grey50'),
          strip.text.x.top = element_text(
            angle = 0, hjust = 0, size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12, hjust = 0.5),
          panel.grid.major.x = element_blank(),
          legend.position = c(0.5, 0.04),
          legend.title.position = 'left',
          legend.background = element_rect(
            fill = "white",        # White background
            color = "black",       # Black border
            linewidth = 0.3        # Border thickness
          )) +
    guides(fill = guide_legend(nrow = 1)) 
}

plot_dist_gap(pw_dist_gap_eval.df)

ggsave('Out/memoire/beta_diff_bray.pdf', 
       bg = 'white', width = 2900, height = 1200, 
       units = 'px', dpi = 200)

ggsave('Out/ISMB2025/beta_diff_bray.pdf', 
       bg = 'white', width = 2500, height = 1000, 
       units = 'px', dpi = 200)

plot_dist_gap(pw_dist_gap_ctrl.df)
ggsave('Out/memoire/beta_diff_bray_ctrls.pdf', 
       bg = 'white', width = 2600, height = 1200, 
       units = 'px', dpi = 200)



## ISMB:
# db_pairs_eval <- list(
#   `A. Kraken 0.45 :\nGTDB 220 – RefSeq` = c('KB45_GTDB', 'KB45'),
#   `B. Sourmash :\nGTDB 220 – RefSeq` = c('SM_gtdb-rs220-rep', 'SM_RefSeq_20250528'),
#   `C. GTDB 220 :\nSourmash – Kraken 0.45` = c('SM_gtdb-rs220-rep', 'KB45_GTDB'),
#   `D. RefSeq :\nSourmash – Kraken 0.45` = c('SM_RefSeq_20250528', 'KB45'),
#   `E. DNA-to-Marker tools :\nmOTUs3 – MetaPhlAn 2023` = c('MOTUS', 'MPA_db2023')
# )



# visualise cv
pairwise_dist_summary %>% 
  filter(!str_detect(Pair_name,"NCBI taxonomy") 
         & Pair_name %in% names(db_pairs_eval)) %>% 
  ggplot(aes(y = Pair_name)) +
  geom_point(size = 3, aes(x = median_diff, fill = Dataset), shape = 23, colour = 'black') +
  geom_point(size = 5, aes(x = mad_diff,colour = Dataset), shape = 4, stroke = 1)

# Are my sets distributed normally?
pw_dist_gap_eval.df %>% 
  group_by(Dataset, Pair_name) %>% 
  slice_sample(n = 5000, replace = FALSE) %>% 
  shapiro_test(dist_diff) %>%
  ggplot(aes(x = p, y = Pair_name, colour = Dataset)) +
  geom_jitter(width = 0, height = 0.2)
# Mostly no!
