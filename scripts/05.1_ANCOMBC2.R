library(pacman)
p_load(magrittr, tidyverse, phyloseq, rlang, vegan, ggridges, ANCOMBC)

ps_species.ls <- read_rds('Out/ps_rare_species.ls.rds')
ps_genus.ls <- read_rds("Out/ps_rare_genus.ls.rds") 
ps_family.ls <- read_rds("Out/ps_rare_family.ls.rds") 

#   
samVar <- 'NAFLD'
ancom_family_NAFLD.ls <- lapply(ps_raw.ls$NAFLD, function(db) {
    ancombc2(
      data = db, 
      tax_level= "Family",
      prv_cut = 0.10, 
      fix_formula = samVar,
      alpha = 0.05,
      verbose = TRUE,
      n_cl = 8)
})

write_rds(ancom_species_NAFLD.ls, 'Out/ancom_species_NAFLD.rds')
write_rds(ancom_genus_NAFLD.ls, 'Out/ancom_genus_NAFLD.rds')
write_rds(ancom_family_NAFLD.ls, 'Out/ancom_family_NAFLD.rds')
write_rds(ancom_genus_RA_Gut.ls, 'Out/ancom_genus_RA_Gut.rds')
write_rds(ancom_species_RA_Gut.ls, 'Out/ancom_species_RA_Gut.rds')
write_rds(ancom_genus_moss.ls,'Out/ancom_genus_moss.rds')
write_rds(ancom_family_moss.ls,'Out/ancom_family_moss.rds')
write_rds(ancom_species_P19_Saliva.ls,'Out/ancom_species_P19_Saliva.ls')
write_rds(ancom_species_P19_Gut.ls,'Out/ancom_species_P19_Gut.ls')

# parse output
compile_ancom <- function(ls, samVar) {
  map(names(ls), function(db) { # iterate over databases
    res <- ls[[db]]$res
    resNames <- names(res)
    sfx <- resNames[str_detect(resNames, paste0("^p_",samVar))] %>% str_remove("^p_")
    res_parsed <- res %>% 
      dplyr::select(taxon, ends_with(sfx)) %>% 
      dplyr::filter(!!sym(paste0("diff_",sfx)) == 1 &
                      !!sym(paste0("passed_ss_",sfx)) == TRUE &
                      !!sym(paste0("q_",sfx)) < 0.05) %>% 
      dplyr::arrange(desc(!!sym(paste0("lfc_",sfx))) ) %>%
      dplyr::mutate(taxon = factor(taxon, levels = unique(taxon))) %>% 
      transmute(LFC = !!sym(paste0("lfc_",sfx)), 
                Taxon = sub(".*?:", "", taxon),
                SE = !!sym(paste0("se_",sfx)),
                adj_p = !!sym(paste0("q_", sfx))) 
    
    res_parsed %>% as_tibble %>% 
      mutate(database = db,
             ref_group = sfx)
    }) %>% list_rbind %>% 
    # Reorder taxa names alphabetically
    mutate(Taxon = fct_relevel(Taxon, sort(levels(Taxon), decreasing = TRUE)))
}

read_rds('Out/ancom_family_NAFLD.rds') %>% 
  compile_ancom('NAFLD') %>% 
  ggplot(aes(x = database, y = Taxon, fill = LFC)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'darkgoldenrod4', mid = 'white', high = 'darkolivegreen3', 
                      na.value = 'grey', midpoint = 0) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )
  
  
  
  
  
  
  