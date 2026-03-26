library(pacman)
p_load(magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
       tidyverse, kableExtra,
       patchwork, ggh4x,
       rstatix)

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
                     'SM_gtdb-rs220-rep', 'SM_RefSeq_20250528', 
                     'MPA_db2023', 'MPA_db2025',
                     'MOTUS3', 'MOTUS4')

# PCoA comparison
permanova.ds <- read_rds('Out/_Rdata/permanova.ds.RDS')%>% 
  filter(Database %in% these_databases) %>% 
  group_by(Dataset) %>%
  mutate(R2_scaled = (R2 - min(R2)) / (max(R2) - min(R2))) %>%
  ungroup() %>% 
  mutate(R2 = round(R2, 3))
 

@@@ formater les décimales à 3
permanova.ds %>% 
  ggplot(aes(x = Database, y = Dataset)) +
  geom_point(aes(size = R2, colour = p)) +
  geom_text(aes(label = R2), colour = 'black', nudge_y = -0.4, size = 3) +
  facet_nested(Dataset ~ CCE_approach + short_alpha_2 + short_alpha_3, scales = 'free') +
  scale_colour_steps2(
    low = "#3d0173", mid = "#f0f0f0", high = "#cccccc",
    midpoint = log10(0.05), 
    trans = "log10",
    breaks = c(0.001, 0.01, 0.05, 0.1, 0.5, 1),
    limits = c(0.001, 1), 
    oob = scales::squish
  ) +
  #scale_size_area(max_size = 10) +
  scale_size(range = c(3, 10)) +
  theme_light()+
  labs(size = "R2", colour = 'p-value') +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = 'grey50'),
    axis.title = element_blank(),
    #panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),   # default is 0.5, shrink gutters between panels
  )

ggsave(paste0('Out/Manuscript/4.beta_permanova.pdf'),
       bg = 'white', width = 2600, height = 1400,
       units = 'px', dpi = 210)


# PLOT variance with p-values
# p_value_lines <- function() {
#   list(
#     geom_vline(
#       aes(xintercept = log10(0.05), linetype = "p = 0.05"), 
#       color = "red", linewidth = 0.2),
#     geom_vline(
#       aes(xintercept = log10(0.01), linetype = "p = 0.01"), 
#       color = "blue", linewidth = 0.2),
#     scale_linetype_manual(
#       name = "p-value thresholds", 
#       values = c("p = 0.05" = "dashed", 
#                  "p = 0.01" = "dashed"))
#   )
# }



# 
# plot_permanova <- function(ds) {
#   ggplot(ds, aes(y = R2, x = log10(p))) +
#     geom_point(
#       shape = 4,  size = 3, stroke = 0.8,
#       aes(colour = Database),
#       position = position_jitter(seed = 1, width = 0.01, height = 0)) +
#     facet_grid(Dataset~., scales = 'free')  +
#     p_value_lines() +
#     scale_colour_manual(values = tooldb_colours, labels = CCE_names) +
#     expand_limits(y = 0) +  
#     scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#     labs(x = expression("log"[10]*"(p-value)"), y = '', colour = 'Methodology') +
#     theme(
#       strip.background = element_rect(fill = 'grey50')
#     )
# }
# 
# p1 <- permanova.ds %>% 
#   filter(
#     Database %in% these_databases,
#     Dataset %in% c('P19_Saliva', 'P19_Gut', 'PD', 'RA_Gut')
#   ) %>% plot_permanova() +
#   labs(y = expression("Proportion of explained variance (R"^2*")"))
# 
# 
# p2 <- permanova.ds %>% 
#   filter(
#     Database %in% these_databases,
#     Dataset %in% c('NAFLD','Bee', 'Moss', 'AD_Skin')
#   ) %>% plot_permanova()
# 
# p1 + p2 + plot_layout(guides = 'collect')
# 
# ggsave('Out/Manuscrit/beta_permanova_bray.pdf', bg = 'white', width = 2200, height = 1200, 
#        units = 'px', dpi = 200)
# 
# # descriptive Statistics
# beta_perm.tbl <- permanova.ds %>% 
#   filter(Database %in% these_databases) %>% 
#   group_by(Dataset) %>% 
#   summarise(
#     minR2 = min(R2),
#     maxR2 = max(R2),
#     n = n()
#   ) %>% 
#   mutate(range = maxR2 - minR2,
#          ratio = maxR2/minR2) %>% 
#   kable(align = "c") %>%
#   kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")); beta_perm.tbl
# 
# save_kable(beta_perm.tbl, 'Out/Manuscrit/tables/beta_perm_table.html')
