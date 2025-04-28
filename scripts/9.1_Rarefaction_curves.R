# Rscript scripts/9.1_Rarefaction_curves.R -c 96 -i Out/ps_full.ls.RDS -o Out/Rarefaction.RDS

library(pacman)
p_load(optparse, tidyverse, phyloseq, rtk, glue, magrittr,
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
  
  make_option(c("-t","--total_cores"), 
              type = "integer", 
              default = 8, 
              help = "Cores to use. If >16, recommend multiples of 16"),
  
  make_option(c("-c","--cores_per_rtk"), 
              type = "integer", 
              default = 2, 
              help = "Cores to use. If >16, recommend multiples of 16")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$input_path)) {
  stop("--input argument is required. See --help for usage.")
}

# Multicore planning
rtk_cores <- opt$cores_per_rtk
list_cores <- floor(opt$total_cores/rtk_cores)
message(glue('Running rtk rarefaction on {list_cores} ps objects in parallel with {rtk_cores} threads each.'))

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
    steps,
    repeats,
    threads) {  # Use one less than available threads
  
  # Extract sequence table
  seqtab <- as(otu_table(ps), 'matrix')
  if(!taxa_are_rows(ps)) { seqtab <- t(seqtab) }
  
  # Remove empty samples (!?)
  sample_depths <- colSums(seqtab) 
  seqtab %<>% .[,which(sample_depths > 0)]
  
  # Calculate depth range
  maxdepth <- max(colSums(seqtab))
  mindepth <- min(colSums(seqtab))
  
  # Depths at which to rarefy
  fractions_of_mindepth <- c(1.1 ,1.2, 1.5, 2, 4, 10, 100, 1000, 2000, 5000, 10000)

  depths <- unique(round(
    c(mindepth/fractions_of_mindepth, # a few below mindepth
      seq(from = mindepth, to = maxdepth, 
          by = floor((maxdepth - mindepth) / (steps - length(fractions_of_mindepth))))
  ))) %>% setdiff(seq(1:10)) # Exclude ultra low depths
  
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
  
  # Extract element lists that correspond to all tested depths 
  rtk_out.ls <- rtk_out.ls[names(rtk_out.ls) %in% rtk_out.ls$depth]
  
  # Extract richness at each depth
  richness_df <- 
    imap(rtk_out.ls, function(rtk_result, depth){
      imap_dfr(rtk_result$divvs, function(sample_out, sample_index) {
        data.frame(
          sample = sample_out$samplename,
          richness = mean(sample_out$richness),
          stringsAsFactors = FALSE,
          depth = as.integer(depth),
          mindepth = mindepth
          )
        })
      }) %>% list_rbind

  # Filter NA results
  return(richness_df)
}

########################
# PARALLEL PROCESSING ###
##########################

# Process from largest to smallest object
all_work <- expand_grid(
  dataset = names(ps.ls),
  database = names(ps.ls[[1]])
  ) %>% mutate(
    size = map2_dbl(dataset, database, ~object.size(ps.ls[[.x]][[.y]]))) %>% 
  arrange(desc(size)) %>% 
  filter(size>0) # Expand grid creates all combinations including non-existing ones

# Process multiple depths in parallel
process_job <- function(row) {
  ps <- ps.ls[[row$dataset]][[row$database]]
  message(glue('Rarefying {row$dataset} + {row$database}...')) ; flush.console() # force messages
  
  rarefaction_out <- rarefaction_curves(
    ps = ps,
    steps = opt$steps,
    repeats = opt$repeats,
    threads = rtk_cores  
  )
  rarefaction_out %>% 
    {
      if (nrow(.) == 0 || ncol(.) == 0) {
        .  # return empty tibble as-is
      } else {
        . %>%
    mutate(Database = row$database, Dataset = row$dataset) %>% 
    filter(!is.na(richness))
      }
    }
  message(glue('Done rarefying {row$dataset} + {row$database}!')); flush.console()
}

# Run with mclapply
results_df <- mclapply(
  split(all_work, 1:nrow(all_work)), 
  process_job,
  mc.cores = list_cores, 
  mc.preschedule = FALSE  # Better for uneven workloads (says deepseek)
) %>% bind_rows() %>% tibble()


# Compute secondary derivatives by sample
result <- results_df %>% 
  group_by(sample) %>% 
  arrange(depth, .by_group = TRUE) %>% 
  mutate(
    first_deriv = (richness - lag(richness)) / (depth - lag(depth)),
    second_deriv = (first_deriv - lag(first_deriv)) / (depth - lag(depth))
  )

# Execute across all list elements 

write_rds(result, opt$output_path)