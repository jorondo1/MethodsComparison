library(pacman)
p_load( magrittr, tidyverse, purrr, phyloseq, patchwork, 
        parallel, vegan)
ps_rare.ls <- read_rds('Out/ps_rare.ls.rds')
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R'))
source("scripts/myFunctions.R")
theme_set(theme_light())

# List of each dataset psmelt (subset of columns)
supermelt <- imap(ps_rare.ls, function(ps_ds.ls, ds) {
  imap(ps_ds.ls, function(ps, db) {
    samvars <- colnames(ps@sam_data)
    psflashmelt(ps) %>% 
      select(Sample, Species, Abundance) %>% 
      group_by(Sample) %>% 
      mutate(relAb = Abundance / sum(Abundance),
             Database = db, .keep = 'unused')
  })# %>% list_rbind 
}) 

example_join <- full_join(
  x = supermelt$P19_Saliva$MPA_db2023,
  y = supermelt$P19_Saliva$MOTUS,
          by = c('Sample', 'Species'))

dim(example_join)
example_join %>%
  filter(!is.na(relAb.x) & !is.na(relAb.y)) %>% 
  group_by(Sample) %>% 
  summarise(comRelAb.x = sum(relAb.x),
            comRelAb.y = sum(relAb.y)) %>% View
