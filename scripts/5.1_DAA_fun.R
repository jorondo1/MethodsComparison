pval_cutoff <- 0.05
ncores <- 4

# Prevalence filtering a phyloseq object
ps_prevalence_filt <- function(ps, threshold) {
  require(phyloseq)
  require(magrittr)
  require(dplyr)
  
  counts.tab <- otu_table(ps)
  
  # Number of samples required to meet threshold
  numSam <- ceiling(threshold * phyloseq::nsamples(ps))
  
  ## Transpose OTU table (species should be arranged by rows)
  if(taxa_are_rows(ps) == FALSE){
    counts.tab <- t(counts.tab)
  }
  
  # Compute taxa prevalence
  keep_taxa <- rowSums(counts.tab != 0) %>% 
    data.frame(count = .) %>% 
    rownames_to_column('taxa') %>% View
    filter(count >= numSam) %$% taxa
  
  res <- phyloseq::prune_taxa(keep_taxa, ps)
  return(res)
}

#############
### ALDEx2 ###
###############
compute_aldex <- function(ps, samVar){
  abund <- ps@otu_table %>% as.data.frame()
  conds <- ps@sam_data %>% as.data.frame() %>% 
    pull(!!sym(samVar)) %>% as.vector()
  
  aldex(reads = abund, 
        conditions = conds,
        mc.samples = 128,
        verbose = TRUE)
}

compile_aldex <- function(results, taxRank, db, ds) {
  results %>% 
    rownames_to_column('Taxon') %>% 
#    dplyr::filter(wi.eBH<pval_cutoff) %>% 
    transmute(Taxon = Taxon, 
           coef = effect,
           adj.p = wi.eBH, 
           taxRank = taxRank,
           database = db,
           dataset = ds,
           DAA_tool = 'Aldex2') %>% as_tibble()
}

##############
### ANCOMBC ###
################
compute_ancombc2 <- function(ps, samVar, taxRank) {
  message(paste('ANCOMBC2 using', ncores, 'cores...'))
 result <- ancombc2(
    data = ps, 
#    tax_level= taxRank,  # DON'T. Since our ps objects have the lowest taxrank as taxa names, we want ANCOM to use that for the output
    prv_cut = 0.10, 
    fix_formula = samVar,
    verbose = TRUE,
    n_cl = ncores)
 result[['samVar']] <- samVar # pass the current samVar to help compilation
 result
}

compile_ancombc2 <- function(results, taxRank, db, ds) {
  resNames <- names(results$res)
  sfx <- resNames[str_detect(resNames, # Extract variable suffix (chosen ref group value)
                             paste0("^p_",results[['samVar']]))] %>% str_remove("^p_")
  
  res_parsed <- results$res %>% 
    dplyr::filter(!!sym(paste0("passed_ss_",sfx)) == TRUE) %>%  # Passed pseudocount sensitivity analysis
    transmute(Taxon = sub(".*?:", "", taxon),
           coef = !!sym(paste0("lfc_",sfx)),
           adj.p = !!sym(paste0("q_", sfx)),
           taxRank = taxRank,
           database = db,
           dataset = ds,
           DAA_tool = 'ANCOMBC2') %>% tibble()
}

##############
### Corncob ###
################
compute_corncob <- function(ps, samVar) {
  my_formula <- as.formula(paste('~', samVar))
  
  corncob::differentialTest(formula= my_formula,
                            formula_null = ~1,
                            phi.formula = ~1,
                            phi.formula_null = ~1,
                            test="Wald", 
                            data=ps,
                            fdr='BH',
                            fdr_cutoff = 0.01,
                            robust = TRUE)
}

compile_corncob <- function(results, taxRank, db, ds) {
  require('purrr'); require('magrittr')
  # Extract coefficients
  coef <- results$all_models %>% 
    purrr::map(~ .x[!is.na(.x)]) %>% 
    compact %>% 
    purrr::map_dfr( \(x) 
                stats::coef(x) %>% 
                  data.frame %>% 
                  rownames_to_column('id') %>% 
                  dplyr::filter(!str_detect(id, '\\(') & 
                                  str_detect(id, 'mu\\.')) %>% 
                  dplyr::select(Estimate)) %$% Estimate
  
  tibble(
    Taxon = names(results$p_fdr),
    coef = coef,
    adj.p = results$p_fdr,
    taxRank = taxRank,
    database = db, 
    dataset = ds,
    DAA_tool = 'corncob'
  ) %>% filter(!is.na(adj.p))
  }

#############
### DESeq2 ###
###############
compute_DESeq2 <- function(ps, samVar) {
  ds_formula <- as.formula(paste('~',samVar))
  dds <- phyloseq_to_deseq2(ps, design = ds_formula)
  dds_res <- DESeq2::DESeq(dds, sfType = 'poscounts')
  results(dds_res, tidy=T, format='DataFrame') %>% tibble()
}

compile_DESeq2 <- function(results, taxRank, db, ds) {
  results %>% 
    # dplyr::filter(!is.na(padj)) %>% # keep significant only
    transmute(Taxon = row,
              coef = log2FoldChange/log2(10),
              adj.p = padj,
              taxRank = taxRank,
              database = db,
              dataset = ds,
              DAA_tool = 'DESeq2')
}

############
### edgeR ###
##############
compute_edgeR <- function(ps, samVar) {
  require('edgeR')
  
  # Compute DAA
  test <- phyloseq_to_edgeR(ps, samVar)
  et = exactTest(test)
  tt = topTags(et, n=nrow(test$table), 
               adjust.method="fdr", 
               sort.by="PValue",
               p.value = 1)
  # Extract results
  tibble(tt@.Data[[1]])
}

compile_edgeR <- function(results, taxRank, db, ds) {
  results %>% 
#    dplyr::filter(FDR < pval_cutoff) %>% # keep significant only
    transmute(Taxon := !!sym(taxRank),
              coef = logFC,
              adj.p = FDR,
              taxRank = taxRank,
              database = db,
              dataset = ds,
              DAA_tool = 'edgeR')
}

###############
### MaAsLin2 ###
#################
# This one is a bit different because the Maaslin2 command creates a written
# output and we will parse from those files. Notice no object need to be created.
compute_Maaslin2 <- function(ps, samVar, taxRank, ds, db, out_path) {
  require('Maaslin2')
  
  maaslin_path <- paste(out_path, '/Maaslin2',taxRank, ds, db, sep = '/')
  
  if (!dir.exists(maaslin_path)) {
    dir.create(maaslin_path, recursive = TRUE)
  }
  
  out <- Maaslin2(
    input_data = otu_table(ps) %>% data.frame(),
    input_metadata = sample_data(ps) %>% data.frame(),
    output = maaslin_path,
    fixed_effects = samVar,
    min_prevalence = 0.1,
    transform = 'LOG',
    standardize = FALSE, 
    normalization = "TSS",
    max_significance = 1,
    plot_heatmap = F, 
    plot_scatter = F,
    cores = detectCores()
  )
}

compile_Maaslin <- function(res_path) {
  Sys.glob(res_path) %>% 
    map_dfr(~ {
      split_string <- str_split(.x, "/", simplify = TRUE)
      read_tsv(.x, col_types = 'cccdddddd') %>%
#        dplyr::filter(qval<pval_cutoff) %>% 
        transmute(Taxon = feature,
                  coef = coef, 
                  adj.p = qval,
                  taxRank = split_string[4],
                  dataset = split_string[5],
                  database = split_string[6],
                  DAA_tool = 'MaAsLin2')
    }) %>% 
    # Maaslin modifies the species names, which is insanely annoying:
    mutate(
      Taxon = str_remove(Taxon, "^\\."),             # 1. Remove leading dot
      Taxon = str_replace(Taxon, "\\.", " "),         # 2. Swap the first dot with a space
      Taxon = str_to_sentence(Taxon),                 # 3. Capitalize the first letter
      Taxon = str_replace_all(Taxon, " sp..", " sp. ")   # 4. Replace '..' with '.'
    ) %>% tibble()
}

############
### RadEMU ###
##############
compute_radEmu <- function(ps, samVar) {
  require('parallel')
  require('radEmu')
  
  # compute fit:
  my_formula <- as.formula(paste('~', samVar))
  ch_fit <- emuFit(formula = my_formula, 
                   Y = ps, run_score_tests = FALSE) 
  ### !DEV! Temporary start**********
  # choose which features have an estimate big enough to have a
  # chance at a significant pval:
  cutoff <- 0
  to_test <- which(abs(ch_fit$coef$estimate) > cutoff)
  ### !DEV Temporary end**********
  
  # Test function for parallelization
  emuTest <- function(category) {
    score_res <- emuFit(
      formula = my_formula,
      Y = ps,
      fitted_model = ch_fit,
      refit = FALSE,
      run_score_tests = TRUE,
      test_kj = data.frame(k = 2, # which cat to test: which("NAFLDNegative" == ch_fit$B %>% rownames) 
                           j = category))
    return(score_res)
  }  # Run in parallel:
  mclapply(to_test, emuTest, mc.cores = ncores)
}

compile_radEmu <- function(results, taxRank, db, ds) {
  map_dfr(results, ~ .x$coef %>%
            dplyr::filter(!is.na(pval)) %>%
            as_tibble()) %>% 
    transmute(
      Taxon = category,
      coef = estimate,
      adj.p = pval,
      taxRank = taxRank,
      database = db,
      dataset = ds,
      DAA_tool = 'radEmu') %>% tibble()  # Add the 'db' name as a new column
}

###############
### ZicoSeq ###
#################
compute_ZicoSeq <- function(ps, samVar) {
  require('rlang')
  require('GUniFrac')
  
  metadata <- ps %>% sample_data() %>% 
    as("data.frame") %>% # convert group var to binary
    mutate(across({{ samVar }}, ~ as.numeric(fct_drop(as.factor(.))) - 1))
  
  result <- ZicoSeq(
    feature.dat = ps %>% otu_table() %>% as.matrix(),
    meta.dat = metadata,
    grp.name = samVar,
    feature.dat.type = 'count',
    #mean.abund.filter = 0.001,
    #prev.filter = 0.1,
    perm.no = 999
  ) 
}

compile_ZicoSeq <- function(results, taxRank, db, ds) {
  coef_matrix <- results$coef.list[[1]]
  tibble(
    Taxon = colnames(coef_matrix),
    coef = coef_matrix[2, ],
    adj.p = results$p.adj.fdr,
    taxRank = taxRank,
    database = db,
    dataset = ds,
    DAA_tool = 'ZicoSeq'
  ) 
}
