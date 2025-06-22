library(pacman)
p_load(tidyverse, slider)
source("scripts/myFunctions.R")

theme_set(theme_light())

rare.df <- readRDS('Out/Rarefaction.rds') %>% 
  filter(depth != 0) %>% 
  # Compute first and second derivatives 
  group_by(Database, sample) %>% 
  arrange(depth, .by_group = TRUE) %>% 
  # Moving average to avoid wild l2f when there's a random slight drop
  
  mutate(mav_richness = slide_dbl(richness, mean, .before = 2, .after = 0, .complete = FALSE)) %>% 
  mutate(
    first_deriv = (richness - lag(richness)) / (depth - lag(depth)),
    second_deriv = (first_deriv - lag(first_deriv)) / (depth - lag(depth)),
    l2f_rate = (mav_richness - lag(mav_richness)) / (log(depth) - lag(log(depth)))
  ) %>% ungroup() %>% 
  mutate(norm_depth = depth/mindepth)
  
# Plot rarefaction curves
these_databases <- c('KB45', 'KB45_GTDB', 'SM_gtdb-rs220-rep', 'SM_RefSeq_20250528')

rare.df %>% 
  #filter(Database %in% c('MOTUS', 'MPA_db2022')) %>% 
  filter(Database %in% these_databases) %>% 
  mutate(Database = factor(Database, levels = these_databases)) %>% 
  # Plot: 
  ggplot(aes(x = depth, y = richness, group = sample, colour = Dataset)) +
  geom_line(linewidth = 0.3, aes(alpha = ifelse(Dataset == "PD", 0.2, 1)))+
  guides(alpha = 'none')+
  #geom_smooth(method = 'loess', se = FALSE, linewidth = 0.1) +
  scale_x_continuous(trans = "log10") +
  facet_grid(Database~., scales = 'free') +
  scale_colour_brewer(palette = 'Set2',
                      guide = guide_legend(override.aes = list(linewidth = 2))) 

# Depth normalized to smallest sample
these_databases <- c('KB10', 'KB45', 'KB45_GTDB', 
                     'SM_gtdb-rs220-rep', 'SM_RefSeq_20250528', 
                     'MPA_db2023','MOTUS')
these_datasets <- c('P19_Gut', 'P19_Saliva', 'PD','NAFLD','AD_Skin', 'Moss', 'Bee') 

rare_plot.df <- rare.df %>% 
  group_by(sample, Database) %>% 
  mutate(z_score = scale(l2f_rate),  
         is_outlier = abs(z_score) > 1
  ) %>% 
  # Remove extreme values (l2f rate is extremely sensitive to 
  # rarefaction randomness, even very small drops in richness cause spikes)
  select(-z_score) %>% 
  ungroup() %>% 
  filter(Database %in% these_databases
         & Dataset %in% these_datasets
         ) %>% 
  mutate(Database = factor(Database, levels = these_databases),
         Dataset = factor(Dataset, levels = these_datasets)) %>% 
   filter(norm_depth < 25)

rare_plot.df %>% 
  ggplot(aes(x = norm_depth, y = richness, colour = Database, 
             group = interaction(sample,Database))) +
  geom_vline(aes(xintercept = 1), color = "black", linewidth = 0.3) +
  geom_line(size = 0.2) +
  # geom_point(size = 2, shape = 8, colour = 'black', alpha = 0.2)+
  facet_grid(Database ~Dataset, scales = 'free') +
  scale_colour_manual(values = tooldb_colours, labels = CCE_names,
                      guide = guide_legend(override.aes = list(linewidth = 2))) +
  labs(x = 'Ratio de profondeur par rapport au plus petit échantillon') +
  theme(legend.position = c(1.04,.5),
        panel.spacing = unit(0.2,'cm'),
        strip.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(r = 3.5, unit = 'cm'),
        legend.background = element_rect(
          fill = "white",        # White background
          color = "black",       # Black border
          linewidth = 0.5        # Border thickness
        )) +
  ylim(0,NA)

ggsave('Out/memoire/rarefaction_curves.pdf',
       bg = 'white', width = 2600, height = 1800,
       units = 'px', dpi = 220)

# Estimate discovery rates at rarefaction depth
rare_plot.df %>%
  filter(!is.na(l2f_rate) & !is_outlier)  %>% 
  # PLOT :
  ggplot(aes(x = norm_depth, y = l2f_rate, colour = Database, group = interaction(Dataset, sample))) +
  #geom_smooth(se = FALSE, linewidth = 0.2) +
  geom_vline(aes(xintercept = 1), color = "black", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0), color = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.2) +
  facet_grid(Database~Dataset, scales = 'free') +
  labs(y = 'Species discovered by doubling depth',
       x = 'Rarefaction depth (ratio to smallest sample in dataset)') +
  scale_colour_manual(values = tooldb_colours, labels = CCE_names,
                      guide = guide_legend(override.aes = list(linewidth = 2))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = c(.92,.3),
        legend.background = element_rect(
          fill = "white",        # White background
          color = "black",       # Black border
          linewidth = 0.5        # Border thickness
        ), 
        strip.text.y = element_blank()) 

ggsave('Out/memoire/rarefaction_l2f_rate.pdf',
       bg = 'white', width = 2600, height = 1800,
       units = 'px', dpi = 180)

# Compare minimum discovery rate with rate at minimum
rare.df %>%
  filter(Dataset == 'P19_Saliva' & Database == 'KB45_GTDB' ) %>% 
  group_by(sample) %>% 
  summarise(max_rate = max(first_deriv, na.rm = TRUE),
            min_rate = min(first_deriv, na.rm = TRUE),
            max_depth = max(depth, na.rm = TRUE)) %>% 
  mutate(lowest_max_depth = min(max_depth)) %>% 
  left_join(
    rare.df %>% 
      select(sample, depth, first_deriv),
    by = "sample"
  ) %>%
  group_by(sample) %>% 
  summarise(
    max_rate = first(max_rate),
    min_rate = first(min_rate),
    # Find the row where depth is closest to lowest_max_depth
    rate_at_min = {
      target_depth = first(lowest_max_depth)
      closest_row = which.min(abs(depth - target_depth))
      first_deriv[closest_row]
    },
    .groups = "drop"
  ) %>% 
  
  ggplot(aes(x = sample)) +
  geom_point(aes(y = min_rate), shape = 1, size = 3 ) +
  geom_point(aes(y = rate_at_min), colour = 'blue') +
  theme(axis.text.x = element_text(angle = 90)) +
  
  labs(y = 'Discovery rate (')


# We need to check whether rarefying at the depth of the shallowest 
# sample, the discovery rate is similar 

plot_df <- rare.df %>% 
  filter(Dataset == 'NAFLD' &
           !is.na(first_deriv) &
           depth <= mindepth)

# This finds the "optimal" scaling using the ratio of largest 
# first to second derivative
axis_scale_factor <- 0.7*( # arbitrary factor for 
  max(abs(plot_df$first_deriv), na.rm = TRUE)/
    max(abs(plot_df$second_deriv), na.rm = TRUE)
) 

# !!! PLOT 2nd derivative convergence
plot_df %>% 
  ggplot(aes(x = log10(depth), colour = sample)) +
  # left y axis: first derivative
  facet_grid(Database~., scales = 'free') +
  geom_line(aes( y = first_deriv)) +
  # right y axis: second derivative (scaled)
  geom_line(data = plot_df,
                          aes( y = axis_scale_factor*second_deriv), linetype = 'dashed') +
  # geom_point(data = filter(plot_df, !is.na(second_deriv) &
  #                          depth == max(depth)), 
  # aes(x = log10(depth), y = second_deriv, colour= sample), 
  # size = 5, shape = 8) +
  theme(legend.position = "none") +
  # Cheat the 2nd y axis back to the original scale
  scale_y_continuous(
    name = "First Derivative (full lines)",
    sec.axis = sec_axis(~ . / axis_scale_factor, name = "Second Derivative (dashed lines)")
  )




