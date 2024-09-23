# intersect in taxonomy
# pairwise Jaccard index on detected species -> similarity matrix
# First test at genus level
library(pacman)
p_load(magrittr, tidyverse, phyloseq)

# 1. Extract taxonomic tables + convert to P/A

ps_genus.ls <- read_rds("Out/ps_rare_genus.ls.rds") 
ps_family.ls <- read_rds("Out/ps_rare_family.ls.rds") 

extract_taxnames <- function(ps, ds, db1, db2, taxRank) {
  
   extract_tax <- function(db) {
     ps[[ds]][[db]] %>% psmelt %>% 
       group_by({{ taxRank }}) %>%
       summarise(Abundance = sum(Abundance)) %>% 
       pull({{ taxRank }}) %>% unique
   }
  tax1 <- extract_tax(db1)
  tax2 <- extract_tax(db2)
  intersect_set <- intersect(tax1,tax2)
  message(paste(length(intersect_set), 'found in common'))
  
  union_set <- union(tax1, tax2)
  message(paste(length(union_set), 'found in total'))
  message(paste(length(intersect_set)/length(union_set), 'intersect'))
  
  setdiff(union_set, intersect_set) %>% sort
}

extract_taxnames(ps_family.ls,
                 'Saliva', 
                 'KB51', 'MOTUS', 
                 Family)

