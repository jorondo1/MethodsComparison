
#############
### ALDEx2 ###
###############
compute_aldex <- function(ps, samVar){
  abund <- ps@otu_table %>% as.data.frame
  conds <- ps@sam_data %>% as.data.frame %>% 
    pull(!!sym(samVar)) %>% as.vector
  
  aldex(reads = abund, 
        conditions = conds,
        mc.samples = 128,
        verbose = TRUE)
}

compile_aldex <- function(results, taxRank, db, ds) {
  results %>% 
    rownames_to_column('Taxon') %>% 
    dplyr::select(effect, wi.eBH, Taxon) %>% 
    dplyr::filter(wi.eBH<0.05) %>% 
    mutate(taxRank = taxRank,
           database = db,
           dataset = ds)
}

##############
### ANCOMBC ###
################

compute_ancombc2 <- function(ps, samVar, taxRank) {
  message(paste('ANCOMBC2 using', ncores, 'cores...'))
 result <- ancombc2(
    data = ps, 
    tax_level= taxRank,
    prv_cut = 0.25, 
    fix_formula = samVar,
    alpha = 0.05,
    verbose = TRUE,
    n_cl = ncores)
 result[['samVar']] <- samVar # pass the current samVar to help compilation
 result
}

compile_ancombc2 <- function(results, taxRank, db, ds) {
  res <- results$res
  samVar <- results[['samVar']]
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
    mutate(taxRank = taxRank,
           database = db,
           dataset = ds)
}

############
### RadEMU ###
##############

compute_radEmu <- function(ps, samVar) {
  require('parallel')
  # compute fit:
  formula_string <- as.formula(paste('~', samVar))
  ch_fit <- emuFit(formula = formula_string, 
                   Y = ps, run_score_tests = FALSE) 
  ### !DEV! Temporary start**********
  # choose which features have an estimate big enough to have a
  # chance at a significant pval:
  cutoff <- 2 
  to_test <- which(abs(ch_fit$coef$estimate) > cutoff)
  ### !DEV Temporary end**********
  
  # Test function for parallelization
  emuTest <- function(category) {
    score_res <- emuFit(
      formula = formula_string,
      Y = ps,
      fitted_model = ch_fit,
      refit = FALSE,
      run_score_tests = TRUE,
      test_kj = data.frame(k = 2, # which cat to test: which("NAFLDNegative" == ch_fit$B %>% rownames) 
                           j = category))
    return(score_res)
  }  # Run in parallel:
  message('Running radEmu...')
  mclapply(to_test, emuTest, mc.cores = ncores)
}

compile_radEmu <- function(results, taxRank, db, ds) {
  map_dfr(results, ~ .x$coef %>%
    dplyr::filter(pval < 0.05) %>%
    as_tibble) %>% 
    transmute(
      coef = estimate,
      pval = pval,
      taxRank = taxRank,
      database = db,
      dataset = ds)  # Add the 'db' name as a new column
}

############
### edgeR ###
##############

compute_edgeR <- function(ps, samVar) {
  # Compute DAA
  test <- phyloseq_to_edgeR(ps, samVar)
  et = exactTest(test)
  tt = topTags(et, n=nrow(test$table), 
               adjust.method="FDR", 
               sort.by="PValue",
               p.value = 0.05)
  # Extract results
  tt@.Data[[1]] %>% tibble
}

compile_edgeR <- function(results, taxRank, db, ds) {
  results %>% 
    dplyr::select(!!sym(taxRank), logFC, FDR) %>% 
    dplyr::filter(FDR < 0.001) %>% # keep significant only
    mutate(database = db, # generic taxa name
           taxRank = taxRank,
           dataset = ds,
           Taxon := !!sym(taxRank), .keep = 'unused') # add database name
}

###############
### MaAsLin2 ###
#################
# This one is a bit different because the Maaslin2 command creates a written
# output and we will parse from those files. Notice no object need to be created.

compute_Maaslin2 <- function(ps, samVar, taxRank, ds, db) {
  maaslin_path <- paste('Out/DAA/Maaslin2',taxRank, ds, db, sep = '/')
  
  if (!dir.exists(maaslin_path)) {
    dir.create(maaslin_path, recursive = TRUE)
  }
  
  out <- Maaslin2(
    input_data = ps %>% otu_table %>% data.frame,
    input_metadata = ps %>% sample_data %>% data.frame,
    output = maaslin_path,
    fixed_effects = samVar,
    transform = 'AST',
    standardize = FALSE, 
    plot_heatmap = F, 
    plot_scatter = F,
    cores = detectCores()-1
  )
}

compile_Maaslin <- function(res_path) {
  Sys.glob(res_path) %>% 
    map_dfr(~ {
      split_string <- str_split(.x, "/", simplify = TRUE)
      read_tsv(.x, show_col_types = FALSE) %>%
        transmute(taxRank = split_string[4],
                  dataset = split_string[5],
                  database = split_string[6],
                  Taxon = feature,
                  coef = coef, 
                  qval = qval)
    }) %>%  dplyr::filter(qval<0.05) %>% 
    # Maaslin modifies the species names, which is insanely annoying:
    mutate(
      Taxon = str_remove(Taxon, "^\\."),             # 1. Remove leading dot
      Taxon = str_replace(Taxon, "\\.", " "),         # 2. Swap the first dot with a space
      Taxon = str_to_sentence(Taxon),                 # 3. Capitalize the first letter
      Taxon = str_replace_all(Taxon, " sp..", " sp. ")   # 4. Replace '..' with '.'
    )
}