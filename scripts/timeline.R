library(pacman)
p_load(tidyverse, lubridate, scales, ggrepel, magrittr)

# Assuming your data is in a TSV with column names : startdate, enddate, side, description, category
# Example structure:
# startdate enddate side  description category
# 2024-01-01  2024-04-30  left  Démonstrateur BIO106 (TPs)  Teaching
# 2024-01-09  2024-01-10  left  Enseigner un atelier Plots with R (Biologie)  Workshop

############################
#======= DATA PREP ========#
############################

timeline_data <- read_tsv("/Users/jorondo/Library/CloudStorage/OneDrive-FreigegebeneBibliotheken–USherbrooke/Isabelle\ Laforest-Lapointe\ -\ RONDEAU_LECLAIRE_Jonathan/Memoire/Comite/timeline.txt") %>%
  mutate(
    description = stringr::str_wrap(description, width = 40),
    
    startdate = as.Date(startdate),
    enddate = as.Date(enddate),
    # Create midpoint for labelling
    midpoint = startdate + (enddate - startdate)/2,
    # Convert to numeric for plotting
    date_numeric = as.numeric(startdate),
    # Create duration column
    duration = as.numeric(enddate - startdate) + 1,
    # Flag for point vs interval
    event_type = ifelse(duration <= 2, "Point", "Interval")
  ) 

# Define your colour categories
category.lvls = c(
  'Milestones' = "#a40000",
  'Conférence' = "#da7901",
  'Formation' = "#007e2f",
  'Enseignement' = "#16317d",
  'Présentations orales' = "#00b7a7",
  'Autre' = "#b86092")

# Reorder factor in the order above:
timeline_data %<>% 
  mutate(category = factor(category,
                           levels = names(category.lvls)))

# Dummy data for unified legend
dummy_legend <- data.frame(
  category = unique(timeline_data$category),
  x = 0, 
  y = as.Date("2024-06-15"),  # Midpoint date that won't be visible
  label = ""  # Empty labels
)

########################
#======= SETUP ========#
########################

bracket_linewidth = 0.2 # thickness of bracket for interval-type events
dash_linewidth = 0.2 # thickness of dash lines
bracket_buffer = 0.1 # how far the bracket tips are from  central line
bracket_height = 0.2 # horizontal height of the brackets
text_position = 0.5 # relative to central line
point_buffer = 0.2 # point-event distance from central line
bracket_xend = bracket_buffer + bracket_height # position of bracket tine root
description_fontsize = 3.2
text_nudge = 0.2 # how far the description text is from the tip of the dashed line.

#######################
#======= PLOT ========#
#######################

ggplot(timeline_data) +
  
  geom_point(
    data = dummy_legend,
    aes(x = x, y = y, color = category),
    size = 0, 
    alpha = 0,
    show.legend = TRUE
  ) +
  
  # Add year markers
  geom_text(data = distinct(timeline_data, year = year(startdate)) %>%
              mutate(year_start = as.Date(paste0(year, "-01-15"))),  # Position in mid-January
            aes(x = 0, 
                y = year_start-22, # fine-tune the year label position vertically
                label = year), 
            size = 3.5, fontface = "bold", 
            color = "gray30") +
  
  # Add month markers
  annotate("text", x = 0, 
           y = seq(min(timeline_data$startdate), 
                   max(timeline_data$enddate)+20, # 20 buffer to add the last months
                   by = "1 month"),
           label = format(seq(min(timeline_data$startdate), 
                              max(timeline_data$enddate)+20, 
                              by = "1 month"), "%b"),
           size = 3, color = "gray40",
           ) +
  
  # === INTERVAL EVENTS - LEFT SIDE
  
  # Lower horizontal segment
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "left"),
               aes(x = -bracket_buffer, xend = -bracket_xend,
                   y = startdate, yend = startdate,
                   color = category),
               size = bracket_linewidth) +
  
  # Upper horizontal segment
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "left"),
               aes(x = -bracket_buffer, xend = -bracket_xend, 
                   y = enddate, yend = enddate,
                   color = category),
               size = bracket_linewidth) +
  
  # Vertical segment
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "left"),
               aes(x = -bracket_xend, xend = -bracket_xend, 
                   y = startdate, yend = enddate,
                   color = category),
               size = bracket_linewidth) +
  
  # Dashed line to event desc
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "left"),
               aes(x = -bracket_xend, xend = -text_position,
                   y = midpoint, yend = midpoint,
                   color = category), 
               size = bracket_linewidth, linetype = "dotted") +
  
  # === INTERVAL EVENTS - LEFT SIDE
  
  # Lower horizontal segment
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "right"),
               aes(x = bracket_buffer, xend = bracket_xend, 
                   y = startdate, yend = startdate,
                   color = category),
               size = bracket_linewidth) +
  
  # Upper horizontal segment
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "right"),
               aes(x = bracket_buffer, xend = bracket_xend, 
                   y = enddate, yend = enddate,
                   color = category),
               size = bracket_linewidth) +
  
  # Vertical segment
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "right"),
               aes(x = bracket_xend, xend = bracket_xend, 
                   y = startdate, yend = enddate,
                   color = category),
               size = bracket_linewidth) +
  
  # Dashed line to event desc
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "right"),
               aes(x = bracket_xend, xend = text_position,
                   y = midpoint, yend = midpoint,
                   color = category),
               size = bracket_linewidth, linetype = "dotted") +
  
  # === POINT EVENTS (both sides)
  
  geom_point(data = filter(timeline_data, event_type == "Point"),
             aes(x = ifelse(side == "left", -point_buffer, point_buffer), 
                 y = startdate, 
                 color = category),
             size = 1, shape = 8) +
  
  geom_segment(data = filter(timeline_data, event_type == "Point"),
               aes(x = ifelse(side == "left", -point_buffer, point_buffer),
                   xend = ifelse(side == "left", -text_position, text_position),
                   y = startdate,
                   yend = startdate,
                   color = category),
               size = dash_linewidth, linetype = "dotted") +
  
  # === DESCRIPTION LABELS + SMART POSITIONING
  
  geom_text_repel(data = filter(timeline_data, side == "left"),
                  aes(x = -text_position, y = midpoint, 
                      label = description, color = category),
                  size = description_fontsize, 
                  hjust = 1, 
                  nudge_x = -text_nudge, 
                  segment.linetype = 'dotted',
                  segment.size = dash_linewidth,
                  max.iter = 100) +
  
  geom_text_repel(data = filter(timeline_data, side == "right"),
                  aes(x = text_position, y = midpoint, 
                      label = description, color = category),
                  size = description_fontsize, 
                  hjust = 0, 
                  nudge_x = text_nudge, 
                  segment.size = dash_linewidth,
                  segment.linetype = 'dotted',
                  max.iter = 100) +
  theme_void() +
  
  # === VISUAL 
  scale_x_continuous(limits = c(-2, 2)) + # enclose a bit
  scale_y_continuous(trans = c('date', "reverse")) + # Reverse y axis
  scale_colour_manual(values = category.lvls, # custom colours
                      labels = names(category.lvls)) +
  theme(
    axis.text = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10)
  ) +
  # Required for tidy legend... remove it and see what happens (:
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 15,        # Squares for all
        size = 6,          # Consistent size
        linetype = "blank", # No lines
        alpha = 1         # No transparency
      )
    )
  ) + labs(color = '') 

ggsave('Out/comite2/timeline.pdf', bg = 'white', width = 1800, height = 2200, 
       units = 'px', dpi = 220)
