library(pacman)
p_load(magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
       tidyverse, kableExtra,
       patchwork,
       rstatix)

ps_rare.ls <- read_rds('Out/_Rdata/ps_rare.ls.rds')
source("scripts/0_Config.R")
theme_set(theme_light())

#   $$\ $$\         $$$$$$$$\ $$$$$$\  $$$$$$\            $$\   $$\ 
#   $$ \$$ \        $$  _____|\_$$  _|$$  __$$\           $$ |  $$ |
# $$$$$$$$$$\       $$ |        $$ |  $$ /  \__|          $$ |  $$ |
# \_$$  $$   |      $$$$$\      $$ |  $$ |$$$$\           $$$$$$$$ |
# $$$$$$$$$$\       $$  __|     $$ |  $$ |\_$$ |          \_____$$ |
# \_$$  $$  _|      $$ |        $$ |  $$ |  $$ |                $$ |
#   $$ |$$ |        $$ |      $$$$$$\ \$$$$$$  |$$\             $$ |
#   \__|\__|        \__|      \______| \______/ \__|            \__|
#

#########################
# PERMANOVA tests ###
###########################

these_databases <- c('KB10', 'KB10_GTDB','KB45', 'KB90', 
                     'KB45_GTDB', 'KB90_GTDB', 
                     'SM_gtdb-rs214-rep', 'SM_gtdb-rs220-rep', 'SM_RefSeq_20250528', 
                     'MPA_db2022','MPA_db2023','MOTUS')

# PCoA comparison
permanova.ds <- read_rds('Out/_Rdata/permanova.ds.RDS')

# PLOT variance with p-values
p_value_lines <- function() {
  list(
    geom_vline(
      aes(xintercept = log10(0.05), linetype = "p = 0.05"), 
      color = "red", linewidth = 0.2),
    geom_vline(
      aes(xintercept = log10(0.01), linetype = "p = 0.01"), 
      color = "blue", linewidth = 0.2),
    scale_linetype_manual(
      name = "p-value thresholds", 
      values = c("p = 0.05" = "dashed", 
                 "p = 0.01" = "dashed"))
  )
}

plot_permanova <- function(ds) {
  ggplot(ds, aes(y = R2, x = log10(p))) +
    geom_point(
      shape = 4,  size = 3, stroke = 0.8,
      aes(colour = Database),
      position = position_jitter(seed = 1, width = 0.01, height = 0)) +
    facet_grid(Dataset~., scales = 'free')  +
    p_value_lines() +
    scale_colour_manual(values = tooldb_colours, labels = CCE_names) +
    expand_limits(y = 0) +  
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = expression("log"[10]*"(p-value)"), y = '', colour = 'Methodology') +
    theme(
      strip.background = element_rect(fill = 'grey50')
    )
}

p1 <- permanova.ds %>% 
  filter(
    Index == 'bray'
    & Database %in% these_databases
    & Dataset %in% c('P19_Saliva', 'P19_Gut', 'PD')
  ) %>% plot_permanova() +
  labs(y = expression("Proportion of inertia explained (R"^2*")"))


p2 <- permanova.ds %>% 
  filter(
    Index == 'bray'
    & Database %in% these_databases
    & Dataset %in% c('NAFLD','Bee', 'Moss', 'AD_Skin')
  ) %>% plot_permanova()

p1 + p2 + plot_layout(guides = 'collect')

ggsave('Out/memoire/beta_permanova_bray.pdf', bg = 'white', width = 2200, height = 1200, 
       units = 'px', dpi = 200)

# descriptive Statistics
beta_perm.tbl <- permanova.ds %>% 
  filter(Index == 'bray' & Database %in% these_databases) %>% 
  group_by(Dataset) %>% 
  summarise(
    minR2 = min(R2),
    maxR2 = max(R2),
    n = n()
  ) %>% 
  mutate(range = maxR2 - minR2,
         ratio = maxR2/minR2) %>% 
  kable(align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")); beta_perm.tbl

save_kable(beta_perm.tbl, 'Out/memoire/tables/beta_perm_table.html')
