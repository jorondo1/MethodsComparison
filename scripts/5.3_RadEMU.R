library(pacman)
p_load(radEmu, purrr)
source('scripts/myFunctions.R')
ps_rare.ls <- read_rds('Out/ps_rare_species.ls.rds')

run_radEmu <- function(ps, taxRank, ds, samVar, ncores = detectCores()-1 ) {
  require('parallel')
  #ps <- ps.ls[[taxRank]][[ds]][[db]]
  # compute fit:
  ch_fit <- emuFit(formula = as.formula(paste('~', samVar)), 
                   Y = ps, run_score_tests = FALSE) 
  # choose which features have an estimate big enough to have a
  # chance at a significant pval:
  cutoff <- 2 # very aggressive for now to facilitate dev
  to_test <- which(abs(ch_fit$coef$estimate) > cutoff)
  # Test function for parallelization
  emuTest <- function(category) {
    score_res <- emuFit(
      formula = ~ NAFLD,
      Y = ps,
      fitted_model = ch_fit,
      refit = FALSE,
      run_score_tests = TRUE,
      test_kj = data.frame(k = 2, # which cat to test: which("NAFLDNegative" == ch_fit$B %>% rownames) 
                           j = category))
    return(score_res)
  }  # Run in parallel:
  score_res <- mclapply(to_test, emuTest, mc.cores = ncores )
}

radEmu_scores.ls <- map(names(ps_rare.ls), function(taxRank) {
  if (taxRank != "Family") {
    return(NULL)  # DEV
  }
  #Process datasets within the "Family" rank
  map(names(ps_rare.ls[[taxRank]]), function(ds) {
    map(names(ps_rare.ls[[taxRank]][[ds]]), function(db) {
      ps <- ps_rare.ls[[taxRank]][[ds]][[db]]
      run_radEmu(ps = ps, taxRank = taxRank,  ds = ds, samVar = group_vars[[ds]])
    })
  })
})
write_rds(radEmu_scores.ls, 'Out/radEmu_scores_Family_NAFLD.rds')

names(radEmu_scores.ls) <- names(ps_rare.ls)
for (taxRank in names(ps_rare.ls)) {
  names(radEmu_scores.ls[['Family']]) <- names(ps_rare.ls[[taxRank]])
}
for (ds in names(ps_rare.ls$Family)) {
  names(radEmu_scores.ls[['Family']][['NAFLD']]) <- names(ps_rare.ls[['Family']][[ds]])
}

compiled_radEmu <- radEmu_scores.ls[['Family']][['NAFLD']] %>%
  imap_dfr(function(db_content, db_name) {
    # Iterate over the elements inside each 'db'
    map_dfr(db_content, ~ .x$coef %>%
              filter(pval < 0.05) %>%
              mutate(database = db_name))  # Add the 'db' name as a new column
  })

# Parse results (weird ass data structure...)
compiled_radEmu %>%
  ggplot(aes(x = database, y = category, fill = estimate)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', mid = 'white', high = 'darkblue', 
                       na.value = 'grey', midpoint = 0) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

# robust_scores <- emuFit(formula = ~ NAFLD, 
#                         Y = ps_rare.ls$Family$NAFLD$MPA_db2023,
#                         fitted_model = ch_fit, refit = FALSE,
#                         run_score_tests = TRUE)

ncores <- detectCores() -1


