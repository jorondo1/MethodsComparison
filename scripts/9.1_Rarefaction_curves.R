library(pacman)
p_load(optparse, tidyverse, phyloseq, rtk, glue, 
       future, furrr, purrr, parallel)

#########################
# +++ ARGUMENT PARSING ###
###########################

option_list <- list(
  make_option(c("-i", "--input_path"), 
              type = "character", 
              default = NULL, 
              help = "Input file path (required)", 
              metavar = "FILE"),
  make_option(c("-o", "--output_path"), 
              type = "character", 
              default = "output_df.rds", 
              help = "Output file path [default: %default]"),
  make_option("--steps", 
              type = "integer", 
              default = 30, 
              help = "Number of rarefaction steps. Will include a step at mindepth and mindepth/c(2,4,6,8,10)"),
  make_option("--repeats", 
              type = "integer", 
              default = 3, 
              help = "Number of times to repeat rarefaction at each step."),
  
  make_option(c("-c","--cores"), 
              type = "integer", 
              default = 2, 
              help = "Cores to use. If >12, recommend multiples of 12")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$input_path)) {
  stop("--input argument is required. See --help for usage.")
}

rtk_cores <- min(12, opt$cores)
list_cores <- floor(opt$cores/rtk_cores)
plan(multisession, workers = list_cores)

message(glue('Running RKT with {rtk_cores} and {list_cores} ps objects in parallel.'))
##################
# +++ LOAD DATA ###
####################

ps.ls <- read_rds(opt$input_path)
ps.ls <- ps.ls$Species

####################
# FUNCTION RAREFY ###
######################

rarefaction_curves <- function(
    ps,
    steps = 50,
    threads = 2,
    repeats = 3) {  # Use one less than available threads
  
  # Extract sequence table
  seqtab <- as(otu_table(ps), 'matrix')
  if(!taxa_are_rows(ps)) { seqtab <- t(seqtab) }
  
  # Calculate depth range
  maxdepth <- max(colSums(seqtab))
  mindepth <- min(colSums(seqtab))
  depths <- round(
    c(mindepth/c(2,4,6,8,10),
      seq(mindepth, maxdepth, floor((maxdepth - mindepth) / (steps - 5)))
  ))
  # Rarefaction, processing samples in parallel
  rtk_out.ls <- suppressWarnings(
    rtk(
      input = seqtab,
      depth = depths,
      repeats = repeats,
      seed = 42,
      threads = threads
      )
  )
  rtk_out.ls <- rtk_out.ls[names(rtk_out.ls) %in% rtk_out.ls$depth]
  
  # Extract richness at each depth
  richness_df <- 
    imap(rtk_out.ls, function(rtk_result, depth){
      imap_dfr(rtk_result$divvs, function(sample_out, sample_index) {
        data.frame(
          sample = sample_out$samplename,
          richness = mean(sample_out$richness),
          stringsAsFactors = FALSE,
          depth = as.integer(depth)
          )
        }) %>% mutate(depth = depth)
      }) %>% list_rbind

  # Filter NA results
  return(richness_df)
}

# Process multiple depths in parallel
results_df <- future_imap(ps.ls, function(ds.ls, dataset) {
  future_imap(ds.ls, function(ps, database) {
    
    message(glue('Rarefying {dataset}, {database}...'))
    
    rarefaction_out <- rarefaction_curves(
      ps = ps,
      steps = opt$steps,
      repeats = opt$repeats,
      threads = rtk_cores
    )
    
    # Add database/dataset columns
    rarefaction_out %>% 
      mutate(
        Database = database,
        Dataset = dataset
      ) %>% tibble() %>% 
      filter(!is.na(richness))
    
  }, .options = furrr_options(seed = TRUE)) %>% list_rbind()
}, .options = furrr_options(seed = TRUE)) %>% list_rbind()

# Reset sequential processing
plan(sequential)

results_df %>% return()
# Combine and format results
result <- results_df %>% 
  tibble() %>% 
  # Compute secondary derivatives by sample
  group_by(sample) %>% 
  arrange(depth, .by_group = TRUE) %>% 
  mutate(
    first_deriv = (richness - lag(richness)) / (depth - lag(depth)),
    second_deriv = (first_deriv - lag(first_deriv)) / (depth - lag(depth))
  )

# Execute across all list elements 

write_rds(result, opt$output_path)