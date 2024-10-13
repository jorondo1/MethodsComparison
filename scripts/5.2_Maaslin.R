library(pacman)
p_load(magrittr, tidyverse, phyloseq, rlang, edgeR, Maaslin2)
source('scripts/myFunctions.R')

ps_rare.ls <- list()
ps_rare.ls[['Species']] <- read_rds('Out/ps_rare_species.ls.rds')
ps_rare.ls[['Genus']]<- read_rds("Out/ps_rare_genus.ls.rds") 
ps_rare.ls[['Family']] <- read_rds("Out/ps_rare_family.ls.rds") 

# function Maaslin
run_Maaslin <- function(ps.ls, taxRank, ds, samVar) {
  require('Maaslin2')
  require('dplyr')
  require('phyloseq')
  require('parallel')
  
  taxRank <- extract_lowest_rank(ps.ls[[taxRank]][[ds]][[1]])
  message(taxRank)
  map(names(ps.ls[[taxRank]][[ds]]), function(db) {
   # dir.create(paste0('Out/Maaslin2',taxRank,,db))
    path <- paste('Out/Maaslin2',taxRank,ds,db, sep = '/')
      
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
    ps <- ps.ls[[taxRank]][[ds]][[db]]
    abund <- data.frame(otu_table(ps))
    metadata <- data.frame(sample_data(ps))
    
    #capture.output(
      out <- Maaslin2(
        input_data = abund,
        input_metadata = metadata,
        output = path,
        fixed_effects = samVar,
        transform = 'AST',
        standardize = FALSE, 
        plot_heatmap = F, 
        plot_scatter = F,
        cores = detectCores()-1
      )
   # )
  })
}

compile_Maaslin <- function(taxRank, ds) {
  require('dplyr')
  require('readr')
  Sys.glob(paste('Out/Maaslin2', taxRank, ds,'*/significant_results.tsv', sep = '/')) %>% 
    map_dfr(~ {
      read_tsv(.x, show_col_types = FALSE) %>%
        transmute(database = basename(dirname(.x)),
                  Taxon = feature,
                  coef = coef, stderr = stderr, qval = qval)
    }) %>%  dplyr::filter(qval<0.05) %>% 
    # Maaslin modifies the species names, which is insanely annoying:
    mutate(
      Taxon = str_remove(Taxon, "^\\."),             # 1. Remove leading dot
      Taxon = str_replace(Taxon, "\\.", " "),         # 2. Swap the first dot with a space
      Taxon = str_to_sentence(Taxon),                 # 3. Capitalize the first letter
      Taxon = str_replace_all(Taxon, " sp..", " sp. ")   # 4. Replace '..' with '.'
    )
}

run_Maaslin(ps_rare.ls, taxRank = 'Family', ds = 'NAFLD', samVar = 'NAFLD')

compile_Maaslin(taxRank = 'Family', ds = 'NAFLD') %>% 
  ggplot(aes(x = database, y = Taxon, fill = coef)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', mid = 'white', high = 'darkblue', 
                       na.value = 'grey', midpoint = 0) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

