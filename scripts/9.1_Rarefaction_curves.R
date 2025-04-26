library(pacman)
p_load(optparse, tidyverse, phyloseq, rtk, parallel, purrr, glue)

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
              help = "Number of times to repeat rarefaction at each step.")
  
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$input_path)) {
  stop("--input argument is required. See --help for usage.")
}

##################
# +++ LOAD DATA ###
####################

ps.ls <- read_rds(input_path)
ps.ls <- ps.ls$Species


####################
# FUNCTION RAREFY ###
######################

rarefaction_curves <- function(
    ps,
    steps = 50,
    repeats = 3, 
    cores = detectCores() - 1) {  # Use one less than available cores
  
  # Extract sequence table
  seqtab <- as(otu_table(ps), 'matrix')
  if(!taxa_are_rows(ps)) { seqtab <- t(seqtab) }
  
  # Calculate depth range
  maxdepth <- max(colSums(seqtab))
  mindepth <- min(colSums(seqtab))
  depths <- c(mindepth/c(2,4,6,8,10),
              seq(mindepth, maxdepth, floor((maxdepth-mindepth)/(steps-5)))
  )
  
  # Parallel processing function
  process_depth <- function(depth, seqtab, repeats) {
    
    # Rarefaction, processing samples in parallel
    rtk_result <- #suppressWarnings(
      rtk(
        input = seqtab,
        depth = depth,
        repeats = repeats,
        seed = 42,
        threads = 2
        )
   # )
    
    # Extract and format results
    richness_df <- tryCatch({
      imap_dfr(rtk_result$divvs, function(sample_out, sample_index) {
        data.frame(
          sample = sample_out$samplename,
          depth = depth,
          richness = mean(sample_out$richness),
          stringsAsFactors = FALSE
        )
      })
    }, error = function(e) {
      data.frame(sample = colnames(seqtab), 
                 depth = depth,
                 richness = NA_real_)
    })
    
    # Filter NA results
    return(richness_df)
  }
  
  # Process multiple depths in parallel
  results <- mclapply(
    depths, 
    process_depth,
    seqtab = seqtab,
    repeats = repeats,
    mc.cores = cores
  )
  
  # Combine and format results
  bind_rows(results) %>% 
    filter(!is.na(richness)) %>% 
    tibble() %>% 
    # Compute secondary derivatives by sample
    group_by(sample) %>% 
    arrange(depth, .by_group = TRUE) %>% 
    mutate(
      depth = log10(depth),
      first_deriv = (richness - lag(richness)) / (depth - lag(depth)),
      second_deriv = (first_deriv - lag(first_deriv)) / (depth - lag(depth))
    ) %>%
    #filter(!is.na(second_deriv)) %>% 
    return()
}

message(glue("Running with {detectCores()} cores!"))

# Execute across all list elements 
results_df <- imap(ps.ls, function(ds.ls, dataset){
  imap(ds.ls, function(ps, database){
    rarefaction_out <- rarefaction_curves(
      ps = ps,
      steps = steps,   
      repeats = repeats,  
      cores = detectCores()     
    )
    
    rarefaction_out %>% 
      mutate(
        Database = database,
        Dataset = dataset
      )
  }) %>% list_rbind
}) %>% list_rbind

write_rds(results_df, output_path)