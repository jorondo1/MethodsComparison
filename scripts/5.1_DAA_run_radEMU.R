# This code was authored by Gemini 2.5 Pro

source('scripts/5.1_DAA_fun_radEMU.R')

# --- Assumed Setup ---
# ps.ls <- your_nested_list
# grouping_variable <- list(PD = "Group", Moss = "Status", ...)

# --- Configuration ---
N_WORKERS <- 40      # Number of parallel jobs to run simultaneously
CORES_PER_JOB <- 2   # Number of cores for each job's mclapply

# --- Execution ---

# 1. Create the flat job table
job_tibble <- create_job_tibble(ps.ls, grouping_variable)
cat("Created a job table with", nrow(job_tibble), "total jobs.\n")
print(head(job_tibble)) # Shows the structure

# 2. Run all jobs in parallel
# The 'func' argument is the name of our newly defined worker function
all_results_df <- run_parallel_jobs(
  job_df = job_tibble,
  func = compute_radEmu_v2,
  n_workers = N_WORKERS,
  cores_per_job = CORES_PER_JOB
)

# 3. Inspect the results
# The output is a tidy tibble with your job specs and a 'results' list-column
print(all_results_df)

# You can now easily filter for failed jobs, for example:
failed_jobs <- all_results_df %>% filter(is.na(results))
print(failed_jobs)

write_rds(all_results_df, 'Out/DAA/radEmu_tmp.RDS')
