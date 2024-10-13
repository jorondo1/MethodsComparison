library(pacman)
p_load(magrittr, tidyverse, phyloseq, rlang, edgeR)

ps_species.ls <- read_rds('Out/ps_rare_species.ls.rds')
ps_genus.ls <- read_rds("Out/ps_rare_genus.ls.rds") 
ps_family.ls <- read_rds("Out/ps_rare_family.ls.rds") 

# function edgeR
compile_edgeR <- function(ps.ls, samVar) {
  require('dplyr')
  require('phyloseq')
  require('edgeR')
  source('https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/phyloseq_to_edgeR.R')

  taxRank <- extract_lowest_rank(ps.ls[[1]])
  # By database: 
  map(names(ps.ls), function(db) {
    # Compute DAA
    test <- phyloseq_to_edgeR(ps.ls[[db]], samVar)
    et = exactTest(test)
    tt = topTags(et, n=nrow(test$table), adjust.method="fdr", sort.by="PValue")
    # Parse results
    tt@.Data[[1]] %>% tibble %>% 
      dplyr::select(!!sym(taxRank), logFC, FDR) %>% 
      dplyr::filter(FDR < 0.05) %>% # keep significant only
      mutate(database = db, # generic taxa name
             Taxon := !!sym(taxRank), .keep = 'unused') # add database name
    }) %>% list_rbind
}

# Plot !
compile_edgeR(ps_family.ls$NAFLD, 'NAFLD') %>% 
  ggplot(aes(x = database, y = Taxon, fill = logFC)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'darkgoldenrod4', mid = 'white', high = 'darkolivegreen3', 
                       na.value = 'grey', midpoint = 0) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )