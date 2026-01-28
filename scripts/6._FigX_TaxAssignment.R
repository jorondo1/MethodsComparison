library(pacman)
p_load( magrittr, mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
        tidyverse, purrr, kableExtra, phyloseq, patchwork, beanplot,
        gghalves,ggh4x,ggridges,
        rstatix, parallel, reshape2, vegan, RColorBrewer, ggridges, htmltools,
        update = FALSE)
ps_rare.ls <- read_rds('Out/_Rdata/ps_rare.ls.rds')
source("scripts/myFunctions.R")
theme_set(theme_light())

################################
# Taxonomic assignment table ####
##################################

# Number of taxa by taxrank at the dataset level
tax_assign_ds <- imap(ps_rare.ls, function(ps_dataset.ls, dataset){
  imap(ps_dataset.ls, function(ps, database){
    tax_count <- ps %>% tax_table() %>% data.frame() %>% tibble()
    
    # Iterate over taxRanks
    map_dfr(taxRanks, function(rank) {
      
      num_tax <- tax_count %>%
        filter(!is.na(!!sym(rank))) %>%
        pull(!!sym(rank)) %>%
        unique() %>%
        length()
      
      tibble(
        Dataset = dataset,
        Database = database,
        Rank = rank,
        Num_tax = num_tax
      ) 
    })
  }) %>% bind_rows() 
}) %>% bind_rows()  %>% 
  mutate(Rank = factor(Rank, levels = taxRanks),
         Database = factor(Database, names(tooldb_colours))) %>% 
  left_join(CCE_metadata, by = 'Database')

write_rds(tax_assign_ds, 'Out/_Rdata/tax_assign_ds.RDS')

# Alt: Number of species by sample
tax_assign_sam <- imap(ps_rare.ls, function(ps_dataset.ls, dataset){
  imap(ps_dataset.ls, function(ps, database){
    
    # Species per sample
    psflashmelt(ps) %>% 
      filter(Abundance > 0) %>% 
      select(Sample, Species) %>% 
      distinct() %>% 
      group_by(Sample) %>% 
      summarise(Num_tax = n()) %>% 
      mutate(
        Dataset = dataset,
        Database = database
      )
  }) %>% bind_rows() 
}) %>% bind_rows()  %>% 
  mutate(Database = factor(Database, names(tooldb_colours))) %>% 
  left_join(CCE_metadata, by = 'Database')

write_rds(tax_assign_sam, 'Out/_Rdata/tax_assign_sam.RDS')

################################
# Taxonomic assignment plot ####
##################################

# Prepare plot data
these_databases <- c('SM_genbank-2022.03', 'SM_gtdb-rs220-rep', 'SM_RefSeq_20250528',
                     'MPA_db2022','MPA_db2023', 'MOTUS',
                     'KB10', 'KB45','KB90', 'KB10_GTDB', 'KB45_GTDB','KB90_GTDB')
these_datasets <- c('Moss', 'NAFLD', 'P19_Gut', 'P19_Saliva', 'PD', 'Bee', 'AD_Skin')

# SUBSET
tax_assign_ds.pdat <- read_rds('Out/_Rdata/tax_assign_ds.RDS') %>% 
  filter(Database %in% these_databases &
           Dataset %in% these_datasets &
           Rank == 'Species') %>% 
  mutate(Prop_db = Num_tax / Num_species_in_db)

tax_assign_sam.pdat <- read_rds('Out/_Rdata/tax_assign_sam.RDS') %>% 
  filter(Database %in% these_databases &
           Dataset %in% these_datasets) %>% 
  mutate(Prop_db = Num_tax / Num_species_in_db)


# Number of Species per datase
tax_assign_ds.pdat %>% 
  ggplot(aes(y = Dataset, x = Num_tax, fill = Database)) +
  geom_col(position = position_dodge2(preserve = 'single',
                                      reverse = TRUE),
           width = 1.1) + # keep bars the same width
  facet_grid(Dataset~., scales = 'free', switch = 'y') +
  scale_fill_manual(values = tooldb_colours, labels = CCE_names) +
  labs(y = 'Dataset', 
       x = "Number of detected species",
       fill = "Methodology") +
  theme(
    panel.spacing.y = unit(0.05, "cm"),
    legend.position = 'inside',
    legend.position.inside = c(.85,.32),
    #axis.text.x = element_text(angle = 45,  hjust=1),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    legend.background = element_rect(
      fill = "white",    
      color = "black",   
      linewidth = 0.5    
    )
  ) 

ggsave('Out/memoire/species_count_ds.pdf',
       bg = 'white', width = 2600, height = 1400,
       units = 'px', dpi = 220)

# Number of Species per sample per dataset
tax_assign_sam.pdat %>% 
  filter(!str_detect(Database, 'KB10')) %>% 
  ggplot(aes(y = "", x = Num_tax, fill = Database)) +
  geom_boxplot(outliers = FALSE, 
               linewidth = 0.3) + 
  facet_wrap(.~Dataset, scales = 'free', ncol = 2) +
  scale_fill_manual(values = tooldb_colours, labels = CCE_names) +
  labs(y = 'Dataset', 
       x = "Number of detected species",
       fill = "Methodology") +
  theme(
    panel.spacing.y = unit(0.05, "cm"),
    legend.position = c(.8,.09),
    #axis.text.x = element_text(angle = 45,  hjust=1),
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    legend.background = element_rect(
      fill = "white",        # White background
      color = "black",       # Black border
      linewidth = 0.5        # Border thickness
    ) 
  )   + guides(fill = guide_legend(ncol = 2)) 


ggsave('Out/memoire/species_count_sam.pdf',
       bg = 'white', width = 2600, height = 1400,
       units = 'px', dpi = 220)


# Min max across all
tax_assign_ds.pdat %>% 
  filter(Taxonomy != 'GTDB') %>% 
  group_by(Dataset) %>% 
  summarise(
    min_species = min(Num_tax),
    max_species = max(Num_tax),
    n = n()
  ) %>% 
  mutate(ratio = max_species / min_species) 

## TABLES 
# Reshape the data for each Dataset
tax_assignment_wide <- tax_assign_ds.pdat %>%
  filter(Database %in% these_databases 
         & Dataset %in% these_datasets
         & Rank == 'Species') %>% 
  select(Dataset, Database, Num_tax) %>% 
  pivot_wider(names_from = Dataset, values_from = Num_tax, values_fill = list(Num_tax = 0))

# Split the data by Dataset and create a kable table for each
tax_assignment_wide %>%
  kable( align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>% 
  save_kable(paste0('Out/memoire/tables/species_count.html'))

# By taxrank, number of taxa found by tool, shown independently by approach
# and within approach, intersect size
# One table per dataset
# Is there a clear break across taxranks ?

## Alternatively : series of Venn diagrams ?
