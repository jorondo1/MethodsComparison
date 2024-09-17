library(pacman)
p_load(phyloseq,tidyverse,cccrm, magrittr,
       car, doParallel, cvequality)

# functions
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/community_functions.R'))
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/rarefy_even_depth2.R'))
source('scripts/myFunctions.R')

objectsToImport <- c("psSalivaKB","psSalivaSM","psSalivaMPA",
                     "psFecesKB","psFecesSM", "psFecesMPA")
for (i in objectsToImport) {assign(i,readRDS(
  paste0("/Users/jorondo/Library/CloudStorage/OneDrive-USherbrooke/Projets/PROVID19/objects/",i,".rds")))}
moss.ps <-readRDS(url('https://github.com/jorondo1/borealMoss/raw/main/data/R_out/mossMAGs.RDS'))

### When doing Feces, remove MSA-3001
ps_list <- list()
meta_parsing <- function(dsName, samData) {
  ps <- list()
  # Metaphlan  
  for (db in c('MPA_db2022', 'MPA_db2023')) {
    ps[[db]] <- parse_MPA(
      MPA_files = paste0(dsName,'/', db, '/*/*profile.txt'),
      column_names = c('Taxonomy', 'NCBI','Abundance', 'Void')) %>% 
      assemble_phyloseq(samData)
  }
  
  # Kraken-bracken (using default headers from parse_MPA function)
  ps[['Bracken51']] <- parse_MPA(
    MPA_files = paste0(dsName,'/Bracken51/*/*_bracken/*_bracken_S.MPA.TXT')) %>% 
    assemble_phyloseq(samData)
  
  # MOTUS
  ps[['MOTUS']] <- parse_MPA(
    MPA_files = paste0(dsName,"/MOTUS/*_profile.txt"), 
    column_names = c('mOTU', 'Taxonomy', 'NCBI', 'Abundance')) %>% 
    assemble_phyloseq(samData)
  # 
  # # Sourmash
  # ps_list[['SM_genbank_202203']] <- left_join(
  #   parse_SM(paste0(dsName,'/SM_genbank_202203/*_genbank-2022.03_gather.csv')),
  #   parse_genbank_lineages(paste0(dsName,'/SM_genbank_202203/genbank-2022.03_lineages.csv')), 
  #   by = 'genome'
  #   ) %>% species_glom() %>% 
  #   assemble_phyloseq(samData)
  # 
  return(ps)
}

ps_list[['Saliva']] <- meta_parsing('Saliva', psSalivaKB@sam_data)
ps_list[['Feces']] <- meta_parsing('Feces', psFecesKB@sam_data)
ps_list[['Moss']] <- meta_parsing('Moss', moss.ps@sam_data)
ps_list$Moss$MPA_db2022 <- NULL
ps_list$Moss$MPA_db2023 <- NULL

## Test with more current MPA database
idx <- c(0,1,2)

# rarefy + diversity 
div.fun <- function(ps) {
  rare.ps <- rarefy_even_depth2(ps, rngseed = 1234) # rarefy
  div_estimate <- list() #initiate list
  for (i in seq_along(idx)) { # compute every diversity index 
    H_q=paste0("H_",i-1) # format H_0, H_1...
    div_estimate[[H_q]] <- estimate_Hill(rare.ps, idx[i])
  }
  return(div_estimate)
}

# Rarefy all datasets and compute diversity
div_rare <- lapply(ps_list, function(sublist) {
  lapply(sublist, div.fun)
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

# Extract tool name and biome name from dataset name
patterns <- str_c(c('saliva'), collapse='|')
Div_long %<>% mutate(
  tool = str_remove(dataset, patterns),
  biome = str_remove(str_remove(dataset, 'ps'), tool),
  # assign abundance estitmation approach
  approach = if_else( # Kraken/Bracken, Sourmash (DNA-to-DNA)
    str_detect(tool, str_c(c('KB', 'SM'), collapse = '|')), "DNA",
    if_else( # Metaphlan, mOTUs (DNA-to-marker)
      str_detect(tool, str_c(c('MPA', 'MOTUS'), collapse = '|')), "marker", 
      NA_character_))
  )

# Visualise differences in diversity across tools
Div_long %>% 
  ggplot(aes(x = database, y = value, fill = database)) +
  geom_boxplot() + geom_line(aes(group = Sample), alpha=0.3)+ theme_light() +
  facet_grid(cols = vars(dataset), rows=vars(index), scales = 'free')

################################################# Not updated below ############
# Pretty obviously yeah
Div_long %>% 
  ggplot(aes(x = tool, y = !!sym(divIDX), fill = tool)) +
  geom_violin() + theme_light()
# model <- lmerTest::lmer(Shannon ~ tool + (1|Sample), data = Shannon_long)
# summary(model)

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
