library(pacman)
p_load(magrittr, tidyverse, phyloseq, rlang, vegan, ggridges, ANCOMBC)

ps_genus.ls <- read_rds("Out/ps_rare_genus.ls.rds") 
ps_family.ls <- read_rds("Out/ps_rare_family.ls.rds") 



### Try to find most variable taxa, because current sample metadata 
# don't yield enough differentially abundant taxa. 

# To leverage all datasets, this function clusters the samples based on the
# top variable taxa (according to clr-transformed abundance) across samples.
clust_sam <- function(ps) {
  clr_wide.mx <- ps %>% 
    microbiome::transform('clr') %>% otu_table
  
  clr_wide <-  clr_wide.mx %>% data.frame %>% 
    rownames_to_column('Taxa')
  
  clr_long <-  clr_wide %>% 
    pivot_longer(names_to = 'Sample', cols = where(is.numeric)) 
  
  clr_taxa_sd <- clr_long %>%  group_by(Taxa) %>% 
    summarise(clr_sd = sd(value)) %>% 
    arrange(desc(clr_sd)) %>% 
    mutate(rank = row_number())
  
  top_taxa <- clr_taxa_sd$Taxa[1:8]
  kmeans_out <- clr_wide %>% 
    filter(Taxa %in% top_taxa) %>% 
    select(-Taxa) %>% t %>% kmeans(centers = 2, iter.max=1000, nstart = 100)
  
  kmeans_out$cluster
}

compile_sample_clust <- function(ls) {
  map(names(ls), function(ps) {
    
    # Generate clusters
    sample_clusters <- clust_sam(ls[[ps]])
    
    # Parse/Compile results
    tibble(
      sample = names(sample_clusters),
      cluster = sample_clusters
    ) %>% 
      mutate(cluster = case_when(cluster == 2 ~ -1,
                                 TRUE ~ cluster))
  }) %>% list_rbind
}

set.seed(34) # 
cluster_assignment <- compile_sample_clust(ps_species.ls$Feces) %>% 
  group_by(sample) %>% 
  summarise(cluster_sum = sum(cluster)) %>% 
  mutate(cluster = case_when(cluster_sum > 0 ~ 'A',
                             TRUE ~ 'B')) %>% 
  column_to_rownames(var = 'sample') 

ps_species_feces.ls <- lapply(ps_species.ls$Feces, function(ps) {
  sample_data(ps) <- cluster_assignment %>% 
    select(-cluster_sum) %>% 
    sample_data %>% 
    cbind(sample_data(ps),.)
  return(ps)
})


#   
samVar <- 'group'
ancom_species_moss.ls <- lapply(ps_species_moss.ls, function(db) {
    ancombc2(
      data = db, 
      tax_level= "Species",
      prv_cut = 0.20, 
      fix_formula = samVar,
      alpha = 0.05,
      verbose = TRUE,
      n_cl = 8)
})
write_rds(ancom_genus_moss.ls,'Out/ancom_genus_moss.rds')
write_rds(ancom_family_moss.ls,'Out/ancom_family_moss.rds')

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
                      !!sym(paste0("q_",sfx)) < 0.01) %>% 
      dplyr::arrange(desc(!!sym(paste0("lfc_",sfx))) ) %>%
      dplyr::mutate(taxon = factor(taxon, levels = unique(taxon))) %>% 
      transmute(LFC = !!sym(paste0("lfc_",sfx)), 
                Taxon = taxon,
                SE = !!sym(paste0("se_",sfx)),
                adj_p = !!sym(paste0("q_", sfx))) 
    
    res_parsed %>% as_tibble %>% 
      mutate(database = db,
             ref_group = sfx)
    }) %>% list_rbind 
}

compile_ancom(ancom_genus_moss.ls, 'Compartment') %>%
  ggplot(aes(x = database, y = Taxon, fill = LFC)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'darkgoldenrod4', mid = 'white', high = 'darkolivegreen3', 
                      na.value = 'grey', midpoint = 0) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

#   dist.mx <- clr_wide.mx %>% t %>% dist
# hc <- hclust(dist.mx)
# clusters <- cutree(hc, k = 2)
  
clr_long %>% filter(Taxa %in% top_taxa) %>% 
  ggplot(aes(x = Sample, y = Taxa, fill = value)) + 
  geom_tile()
  
  
  
  
  
  
  
  