# # This code was authored by Gemini 2.5 Pro
# /// Prompt: 
# here's a function I pass a R list to, as well as a function to be applied.
# 
# compute_3_lvl <- function(ps.ls, func, ...){
#   require('furrr') 
#  ... <function in scripts/myFunctions.R>
# }
# 
# The function to be applied is the fololwing:
#   compute_radEmu <- function(ps, samVar) {
#  ... <function in scripts/myFunctions.R>
#   }
# so the call is 
# test_radEmu <- compute_3_lvl(ps.ls, func = compute_radEmu)
# # 
# Now this function runs a very long time on a single core (or 2 i'm not sure). Usually I use 
# future::plan(multisession, workers = 6)
# ncores=12
# (we have about 90 cores available)
# 
# However I'd like to run way more workers with only 2 cores each. The thing is, as you can see from the 3 level function, the physeq objects I'm using within the list are nested at the 3rd level of the list, e.g. ps.ls$Species$PD$MPA. My code structure can only paralellise the last-level of the imap loop. Since there are multiple sublists at each level, e.g. 
# ps.ls$Genus$Moss$MPA. the current structure doesn't allow running multiple datasets in parallel outside of those nested in the last level of the list (e.g. all ps objects in $PD can be run in parallel, but not at the same time as the objects in $Moss). Maybe I'm wrong. But I'd like to be able to run all objects in parallel at once to saturate the workers (I have about 12 objects per last level, 8 dataset levels (e.g. Moss, PD) and two first level (Species, Genus). I could run with 36 workers at 2 cores each fdor example"
# What would you do? We could either change the ps.ls structure, but I'm using the level names with imap to pass to the final output. Or we can create a new compute_3_lvl function...
# \\\prompt

library(purrr)
library(dplyr)

# Flatten the ps.ls list

#' Create a job specification tibble from the nested phyloseq list.
#'
#' @param ps.ls The 3-level nested list of phyloseq objects.
#' @param grouping_variable A named list mapping dataset names ('ds') to their grouping variable.
#' @return A tibble where each row is a unique job to be run.
create_job_tibble <- function(ps.ls, grouping_variable) {
  
  # Use purrr's mapping functions to walk the nested list and extract jobs
  jobs <- imap_dfr(ps.ls, function(taxRank.ls, taxRank) {
    imap_dfr(taxRank.ls, function(ds.ls, ds) {
      imap_dfr(ds.ls, function(db.ps, db) {
        
        # Create a one-row tibble for this specific combination
        tibble(
          taxRank = taxRank,
          ds = ds,
          db = db,
          ps = list(db.ps),  # Keep the phyloseq object in a list-column
          samVar = grouping_variable[[ds]]
        )
      })
    })
  })
  
  return(jobs)
}


# Have radEmu accept core_per_job
compute_radEmu_v2 <- function(ps, samVar, taxRank, ds, db, cores_per_job = 2) {
  
  # You can use taxRank, ds, db here for logging or custom file paths
  message(sprintf("STARTING JOB: TaxRank=%s, Dataset=%s, DB=%s on %d cores", 
                  taxRank, ds, db, cores_per_job))
  
  # compute fit:
  my_formula <- as.formula(paste('~', samVar))
  
  # It's good practice to wrap in a tryCatch to prevent one failure from
  # stopping the entire parallel run.
  result <- tryCatch({
    ch_fit <- radEmu::emuFit(formula = my_formula, 
                             Y = ps, run_score_tests = FALSE) 
    
    cutoff <- 0
    to_test <- which(abs(ch_fit$coef$estimate) > cutoff)
    
    if (length(to_test) == 0) {
      message("No features to test. Skipping parallel score tests.")
      return(list(fit = ch_fit, scores = NULL))
    }
    
    # Test function for parallelization
    emuTest <- function(category) {
      # This function will run inside mclapply, so it needs radEmu
      radEmu::emuFit(
        formula = my_formula,
        Y = ps,
        fitted_model = ch_fit,
        refit = FALSE,
        run_score_tests = TRUE,
        test_kj = data.frame(k = 2, j = category)
      )
    }
    
    # Run in parallel, using the provided number of cores
    score_results <- parallel::mclapply(to_test, emuTest, mc.cores = cores_per_job)
    
    message(sprintf("COMPLETED JOB: TaxRank=%s, Dataset=%s, DB=%s", taxRank, ds, db))
    
    list(fit = ch_fit, scores = score_results)
    
  }, error = function(e) {
    message(sprintf("ERROR in JOB: TaxRank=%s, Dataset=%s, DB=%s. Error: %s", 
                    taxRank, ds, db, e$message))
    return(NA) # Return NA or an error object on failure
  })
  
  return(result)
}

library(furrr)
library(future)

#' Run a function in parallel over a flattened job tibble.
#'
#' @param job_df A tibble created by create_job_tibble().
#' @param func The function to apply to each row (e.g., compute_radEmu_v2).
#' @param n_workers The number of parallel workers (jobs to run at once).
#' @param cores_per_job The number of cores EACH worker will use for its internal mclapply.
#' @param ... Additional fixed arguments to pass to `func`.
#' @return A list containing the results of each job.

run_parallel_jobs <- function(job_df, func, n_workers, cores_per_job, ...) {
  
  # IMPORTANT: Setup the parallel backend.
  # 'multicore' is best here because your inner task uses mclapply (forking).
  # Total cores used = n_workers * cores_per_job.
  # Make sure this doesn't exceed your machine's capacity.
  # e.g., 36 workers * 2 cores/job = 72 cores.
  plan(multicore, workers = n_workers)
  
  message(sprintf("Starting parallel execution with %d workers, each using %d cores.",
                  n_workers, cores_per_job))
  
  # We use future_pmap to iterate over the rows of the tibble.
  # pmap passes the columns of job_df as named arguments to the function.
  # We also need to pass 'cores_per_job' and any other arguments.
  
  all_results <- future_pmap(
    .l = job_df, 
    .f = func,
    cores_per_job = cores_per_job,
    ...,
    .options = furrr_options(
      seed = TRUE, 
      # Packages that need to be loaded on each forked worker
      packages = c("radEmu", "parallel") 
    ),
    .progress = TRUE
  )
  
  # Add the job specifications to the results for easy identification
  results_with_specs <- bind_cols(job_df %>% select(-ps), tibble(results = all_results))
  
  return(results_with_specs)
}

