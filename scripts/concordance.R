library(pacman)
p_load(phyloseq,tidyverse,cccrm, magrittr,
       car, doParallel, cvequality)

# functions
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/community_functions.R'))
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/rarefy_even_depth2.R'))
source('scripts/myFunctions.R')

# Import External data
objectsToImport <- c("psSalivaKB", "psFecesKB")
for (i in objectsToImport) {assign(i,readRDS(
  paste0("/Users/jorondo/Library/CloudStorage/OneDrive-USherbrooke/Projets/PROVID19/objects/",i,".rds")))}
moss.ps <-readRDS(url('https://github.com/jorondo1/borealMoss/raw/main/data/R_out/mossMAGs.RDS'))

# Function to call parsing functions for each dataset,
# results in a list containing 1 phyloseq object by dataset by tool (database)
# with general structure List$Dataset$Database.ps
meta_parsing <- function(dsName, samData, filtering = FALSE) {
  ps <- list()
  # Metaphlan  
  for (db in c('MPA_db2022', 'MPA_db2023')) {
    ps[[db]] <- parse_MPA(
      MPA_files = paste0(dsName,'/', db, '/*/*profile.txt'),
      column_names = c('Taxonomy', 'NCBI','Abundance', 'Void')) %>% 
      assemble_phyloseq(samData, filtering = filtering)
  }
  
  # Kraken-bracken (using default headers from parse_MPA function)
  ps[['Bracken51']] <- parse_MPA(
    MPA_files = paste0(dsName,'/Bracken51/*/*_bracken/*_bracken_S.MPA.TXT')) %>% 
    assemble_phyloseq(samData, filtering = filtering)
  
    # MOTUS
  ps[['MOTUS']] <- parse_MPA(
    MPA_files = paste0(dsName,"/MOTUS/*_profile.txt"), 
    column_names = c('mOTU', 'Taxonomy', 'NCBI', 'Abundance')) %>% 
    assemble_phyloseq(samData, filtering = filtering)
  
  # Sourmash
  ps[['SM_genbank_202203']] <- left_join(
    parse_SM(paste0(dsName,'/SM_genbank_202203/*_genbank-2022.03_gather.csv')),
    parse_genbank_lineages(paste0(dsName,'/SM_genbank_202203/genbank-2022.03_lineages.csv')),
    by = 'genome'
    ) %>% species_glom() %>%
    assemble_phyloseq(samData, filtering = filtering)
  
  for (db in c(#'SM_gtdb_rs214',
               'SM_gtdb_rs214_rep')) {
    ps[[db]] <- left_join(
      parse_SM(paste0(dsName,'/', db, '/*_gather.csv')),
      parse_GTDB_lineages(paste0(dsName,'/', db, '/', db, '_lineages.csv')),
      by = 'genome'
    ) %>% assemble_phyloseq(samData, filtering = filtering)
  }
  return(ps)
}

ps_list <- list()
ps_list[['Saliva']] <- meta_parsing('Saliva', psSalivaKB@sam_data)
ps_list[['Feces']] <- meta_parsing('Feces', psFecesKB@sam_data)
ps_list[['Moss']] <- meta_parsing('Moss', moss.ps@sam_data)
ps_list[['Moss']][['SM_gtdb_rs214_rep_MAGs']] <- moss.ps 
ps_list$Moss$MPA_db2022 <- NULL
ps_list$Moss$MPA_db2023 <- NULL

ps_list_filtered <- list()
ps_list_filtered[['Saliva']] <- meta_parsing('Saliva', psSalivaKB@sam_data, filtering = TRUE)
ps_list_filtered[['Feces']] <- meta_parsing('Feces', psFecesKB@sam_data, filtering = TRUE)
ps_list_filtered[['Moss']] <- meta_parsing('Moss', moss.ps@sam_data, filtering = TRUE)
ps_list_filtered[['Moss']][['SM_gtdb_rs214_rep_MAGs']] <- moss.ps%>% filter_low_prevalence()
ps_list_filtered$Moss$MPA_db2022 <- NULL
ps_list_filtered$Moss$MPA_db2023 <- NULL

# compute sparseness for all datasets
sparseness.df <- compile_sparseness(ps_list)
sparseness_filtered.df <- compile_sparseness(ps_list_filtered)

# Vizualise :
bind_rows(sparseness.df, sparseness_filtered.df, .id = 'filtered') %>% 
  mutate(filtered = case_when(filtered ==1 ~ 'no', TRUE ~'yes')) %>% 
ggplot(aes(y = sparseness, x = dataset, fill = filtered)) +
  geom_col(position = 'dodge') + theme_minimal() +
  facet_grid(cols=vars(database), scales = 'free')

# Rarefy all datasets
ps_list_rare <- lapply(ps_list_filtered, function(sublist) {
  lapply(sublist, rarefy_even_depth2, rngseed = 1234)
})

# Compute diversity across indices
div_rare <- lapply(ps_list_rare, function(sublist) {
  lapply(sublist, div.fun, idx = c(0,1,2))
})

# Check distribution of indices :
lapply(div_rare, function(sublist) {
  lapply(sublist, function(subsublist) {
    lapply(subsublist, function(element) {
      shapiro.test(element)$p.value
    })})
}) 

# Create long dataframe 
Div_long <- map(names(div_rare), function(ds) { #iterate over dataset names
  ds_sublist <- div_rare[[ds]] 
  map(names(ds_sublist), function(db) { # iterate over databases
    db_sublist <- ds_sublist[[db]]
    map(names(db_sublist), function(hill) { # iterate over index types
      values <- db_sublist[[hill]]
      tibble( # build dataset
        Sample = names(values),
        dataset = ds,
        database = db,
        index = hill,
        value = values
      )
    }) %>% list_rbind # collapse list into single df
  }) %>% list_rbind  
}) %>% list_rbind

# Visualise differences in diversity across tools
Div_long %>% 
  ggplot(aes(x = database, y = value, fill = database)) +
  geom_boxplot() + geom_line(aes(group = Sample), alpha=0.3)+ theme_light() +
  facet_grid(cols = vars(dataset), rows=vars(index), scales = 'free')







################################################# Not updated below ############
# Pretty obviously yeah

# Dynamic formula : 
formula <- reformulate('tool', response = divIDX)
leveneTest(formula, data = Div_long) # same variance
kruskal.test(formula, data = Div_long) # different distributions

#### Coefficient of variation ; testing whether they are different across categories
# Coefficient of variation for each sample: 
apply(Div, 2, function(x) (sd(x)/mean(x))) # Are they really different ?
# Using Feltz and Miller 1996 implemented in : 
# https://cran.r-project.org/web/packages/cvequality/vignettes/how_to_test_CVs.html
mslr_test(1e5, Div_long[[divIDX]], Div_long$tool) # Same variation coefficient

#############################
### Bland & Altman Analysis ###
#############################

BA_analysis <- compute_meandiff(Div, 'KB', 'SM')
(test.norm <- shapiro.test(BA_analysis$Diff) %$% p.value)

# If differences are not normal, log transform diversity counts :
if(test.norm <0.05) {
  BA_analysis <- Div %>% 
    mutate(across(everything(), log)) %>% 
    compute_meandiff('KB', 'SM')
  shapiro.test(BA_analysis$Diff) %$% p.value
}

# global statistics of means and differences:
mean_diff <- mean(BA_analysis$Diff)
s <- sd(BA_analysis$Diff)
upLimit = mean_diff + 1.96*s
loLimit = mean_diff - 1.96*s
n = length(BA_analysis$Diff)
mean_sd = sqrt(s^2/n)
s_sd = sqrt((3*s^2)/n)

# proportion of differences outside limits
mean(BA_analysis$Diff < loLimit | BA_analysis$Diff > upLimit) 

# Bland-Altman Plot :
ggplot(BA_analysis) +
  theme_light() +
  annotate(geom = 'rect', xmin = -Inf, xmax = Inf, ymin = upLimit - s_sd, 
           ymax =  upLimit + s_sd, alpha = 0.2, fill = 'blue') +
  annotate(geom = 'rect', xmin = -Inf, xmax = Inf, ymin = mean_diff - mean_sd, 
           ymax = mean_diff + mean_sd, alpha = 0.2, fill = 'red') +
  annotate(geom = 'rect', xmin = -Inf, xmax = Inf, ymin = loLimit - s_sd, 
           ymax = loLimit + s_sd, alpha = 0.2, fill = 'blue') +
  geom_hline(yintercept = mean_diff) +
  geom_hline(yintercept = upLimit, linetype = 'dashed') +
  geom_hline(yintercept = loLimit, linetype = 'dashed') +
  geom_point(aes(x = mean, y = Diff)) +
  labs(x = "Mean", y = "Difference", title = "Limits of Agreement between Kraken and Sourmash")

#### Trying out something... plot 3-way variation against mean
sd_mean <- Div %>% rowwise %>% 
  transmute(div_mean = mean(c_across(everything())),
            div_sd = sd(c_across(everything()))) 

cor.test(sd_mean$div_mean, sd_mean$div_sd, tool = 'spearman') # no correlation
plot(div_sd~div_mean, data = sd_mean)

####################################################
### (Overall) Concordance Correlation Coefficient ###
# Should we scale the diversities for the CCC test ? Only center?
####################################################

# epiR tool https://search.r-project.org/CRAN/refmans/epiR/html/epi.occc.html
Div %>% select(KB, SM) %>% 
  # rationale for scaling ?
  mutate(across(everything(), ~scale(.x, scale = FALSE))) %>% 
  epiR::epi.occc(pairs = TRUE)

# CCCRM tool
Div_scaled_long <- Div %>% 
  # Center only, because variances are equal
  mutate(across(everything(), scale)) %>% 
  rownames_to_column("Sample") %>% 
  pivot_longer(values_to = divIDX, 
               names_to = "tool",
               cols = all_of(tools)) %>% 
  mutate(tool = factor(tool, levels = tools)) %>% 
  left_join(psSalivaKB@sam_data %>% data.frame %>% 
              dplyr::select(treatDay) %>% 
              rownames_to_column("Sample"),
            by = "Sample")

ccc_result <- cccvc(Div_scaled_long %>% filter(tool !="MPA"), ry = divIDX, 
                    rind = "Sample", rmet = "tool")
summary(ccc_result)

ccc_result <- cccvc(Div_scaled_long %>% 
                      filter(tool != 'KB'), 
                    ry = "Shannon", 
                    rind = "Sample", rmet = "tool")
summary(ccc_result)
