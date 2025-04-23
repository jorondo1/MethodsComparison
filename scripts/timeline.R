library(pacman)
p_load(ggplot2, dplyr, lubridate, scales,ggrepel)

# Assuming your data is in a CSV with columns: startdate, enddate, description, category
# Example structure:
# startdate,enddate,description,category
# 2024-01-01,2024-04-30,Démonstrateur BIO106 (TPs),Teaching
# 2024-01-09,2024-01-10,Enseigner un atelier Plots with R (Biologie),Workshop

# Read data
timeline_data <- read_tsv("/Users/jorondo/Library/CloudStorage/OneDrive-FreigegebeneBibliotheken–USherbrooke/Isabelle\ Laforest-Lapointe\ -\ RONDEAU_LECLAIRE_Jonathan/Memoire/Comite/timeline.txt") %>%
  mutate(
    description = stringr::str_wrap(description, width = 40),
    
    startdate = as.Date(startdate),
    enddate = as.Date(enddate),
    # Create midpoint for labeling
    midpoint = startdate + (enddate - startdate)/2,
    # Convert to numeric for plotting
    date_numeric = as.numeric(startdate),
    # Create duration column
    duration = as.numeric(enddate - startdate) + 1,
    # Flag for point vs interval
    event_type = ifelse(duration <= 2, "Point", "Interval")
  ) 

category.lvls = c(
  'Milestones' = "#a40000",
  'Conférence' = "#da7901",
  'Formation' = "#007e2f",
  'Enseignement' = "#16317d",
  'Présentations orales' = "#00b7a7",
  'Autre' = "#b86092")

timeline_data %<>% 
  mutate(category = factor(category,
                           levels = names(category.lvls)))

bracketsize = 0.2

# Dummy data for unified legend
dummy_legend <- data.frame(
  category = unique(timeline_data$category),
  x = 0, 
  y = as.Date("2024-06-15"),  # Midpoint date that won't be visible
  label = ""  # Empty labels
)


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
            aes(x = 0, y = year_start-24, label = year),
            size = 3.5, fontface = "bold", 
            color = "gray30") +
  
  # Add month markers
  annotate("text", x = 0, 
           y = seq(min(timeline_data$startdate), max(timeline_data$enddate)+20, by = "1 month"),
           label = format(seq(min(timeline_data$startdate), max(timeline_data$enddate)+20, by = "1 month"), "%b"),
           size = 2.8, color = "gray40",
           ) +
  
  
  # Interval events - left side
  ## Lower horizontal segment
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "left"),
               aes(x = -0.1, xend = -0.3, 
                   y = startdate, yend = startdate,
                   color = category),
               size = bracketsize) +
  ## Upper horizontal segment
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "left"),
               aes(x = -0.1, xend = -0.3, 
                   y = enddate, yend = enddate,
                   color = category),
               size = bracketsize) +
  ## Vertical segment
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "left"),
               aes(x = -0.3, xend = -0.3, 
                   y = startdate, yend = enddate,
                   color = category),
               size = bracketsize) +
  ## Dashed line to event desc
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "left"),
               aes(x = -0.3, xend = -0.5,
                   y = midpoint, yend = midpoint,
                   color = category),
               size = bracketsize, linetype = "dotted") +
  
  # Interval events - right side
  ## Lower horizontal segment
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "right"),
               aes(x = 0.1, xend = 0.3, 
                   y = startdate, yend = startdate,
                   color = category),
               size = bracketsize) +
  ## Upper horizontal segment
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "right"),
               aes(x = 0.1, xend = 0.3, 
                   y = enddate, yend = enddate,
                   color = category),
               size = bracketsize) +
  ## Vertical segment
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "right"),
               aes(x = 0.3, xend = 0.3, 
                   y = startdate, yend = enddate,
                   color = category),
               size = bracketsize) +
  ## Dashed line to event desc
  geom_segment(data = filter(timeline_data, event_type == "Interval", side == "right"),
               aes(x = 0.3, xend = 0.5,
                   y = midpoint, yend = midpoint,
                   color = category),
               size = bracketsize, linetype = "dotted") +
  
  # Point events - both sides (modified version)
  geom_point(data = filter(timeline_data, event_type == "Point"),
             aes(x = ifelse(side == "left", -0.2, 0.2), 
                 y = startdate, 
                 color = category),
             size = 1, shape = 8) +
  
  geom_segment(data = filter(timeline_data, event_type == "Point"),
               aes(x = ifelse(side == "left", -0.2, 0.2),
                   xend = ifelse(side == "left", -0.5, 0.5),
                   y = startdate,
                   yend = startdate,
                   color = category),
               size = 0.2, linetype = "dotted") +
  
  # Labels with smart positioning
  geom_text_repel(data = filter(timeline_data, side == "left"),
                  aes(x = -0.5, y = midpoint, 
                      label = description, color = category),
                  size = 3, 
                  hjust = 1, 
                  direction = "y",
                  nudge_x = -0.1, 
                  segment.linetype = 'dotted',
                  segment.size = 0.2,
                  max.iter = 1000) +
  
  geom_text_repel(data = filter(timeline_data, side == "right"),
                  aes(x = 0.5, y = midpoint, 
                      label = description, color = category),
                  size = 3, 
                  hjust = 0, 
                  direction = "y",
                  nudge_x = 0.1, 
                  segment.size = 0.2,
                  segment.linetype = 'dotted',
                  max.iter = 1000) +
  theme_minimal() +
  
  # Visual styling
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_date(breaks = date_breaks("1 month"),
               labels = function(x) format(x, "%b"),  # Just 3-letter month abbreviations
               expand = expansion(add = c(30, 30))) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text = element_text(size = 10, angle = 90, vjust = 0.5),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.spacing.x = unit(5, "pt"),
    legend.text = element_text(size = 9)
  ) +
  scale_colour_manual(values = category.lvls, labels = names(category.lvls)) +
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

ggsave('Out/comite2/timeline.png', bg = 'white', width = 1600, height = 2200, 
       units = 'px', dpi = 200)
