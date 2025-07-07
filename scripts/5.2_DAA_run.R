# module load StdEnv/2023 r/4.4.0
library(pacman)
p_load(magrittr, tidyverse, phyloseq,
       furrr, purrr, parallel,
       ANCOMBC, ALDEx2, corncob, radEmu, 
       edgeR, DESeq2, Maaslin2, GUniFrac,
       update = FALSE)

source('scripts/myFunctions.R')
source('scripts/5.1_DAA_fun.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/tax_glom2.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/rarefy_even_depth2.R')

# ps.ls.tmp <- read_rds('Out/_Rdata/ps_filt.ls.RDS') 

# ps.ls <- list()
# #ps.ls[['Species']] <- ps.ls.tmp
# for (taxRank in c('Genus', 'Family')) {
#   ps.ls[[taxRank]] <- lapply(ps.ls.tmp, function(ds) {
#     lapply(ds, function(db) {
#       tax_glom2(db, taxrank = taxRank)
#     })
#   })
# }
# 
# # Rarefied dataset for MetaPhlAn
# ps_rare.ls <- lapply(ps.ls, function(taxRank) {
#   lapply(taxRank, function(ds) {
#     lapply(ds, function(db) {
#       rarefy_even_depth2(db, rngseed = 42, ncores = 24)
#     })
#   })
# })
# 
# write_rds(ps.ls, 'Out/_Rdata/ps_taxranks.ls.RDS')
# write_rds(ps_rare.ls, 'Out/_Rdata/ps_rare_taxranks.ls.RDS')

remove_databases <- function(ps.list, to_remove) {
  lapply(ps.list, function(taxRank) {
    lapply(taxRank, function(ds){
      imap(ds, function(dataset, db) {
        if (db %in% to_remove) {
          NULL  # Replace with NULL if name matches
        } else {
          dataset  # Keep if no match
        }
      }) %>% compact() #remove nulls 
    })
  }) 
}

ps.ls <- remove_databases(
  read_rds('Out/_Rdata/ps_taxranks.ls.RDS'),
  c("SM_genbank-2022.03", "KB10", "KB10_GTDB"))
ps_rare.ls <- remove_databases(
  read_rds('Out/_Rdata/ps_rare_taxranks.ls.RDS'),
  c("SM_genbank-2022.03", "KB10", "KB10_GTDB"))

# To work on a subset at whichever list level: 
# ps.ls <- map(ps.ls, ~ .x["NAFLD"])
# ps_rare.ls <- map(ps_rare.ls, ~ .x["NAFLD"])

out_path <- 'Out/DAA'
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

future::plan(multisession, workers = 14)
ncores=6
# Aldex2
test_aldex <- compute_3_lvl(ps.ls, func = compute_aldex)
compile_3_lvl(test_aldex, func = compile_aldex) %>% 
  write_tsv(paste0(out_path,'/Aldex2.tsv')) 

# ANCOM-BC2
test_ancombc2 <- compute_3_lvl(ps.ls, func = compute_ancombc2)
compile_3_lvl(test_ancombc2, func = compile_ancombc2) %>% 
  write_tsv(paste0(out_path,'/AncomBC2.tsv'))

# Corncob
test_corncob <- compute_3_lvl(ps.ls, func = compute_corncob)
compile_3_lvl(test_corncob, func = compile_corncob) %>% 
  write_tsv(paste0(out_path, '/corncob.tsv'))

# DESeq2 !Note: L2FC transformed to L10FC for scale
test_DESeq2 <- compute_3_lvl(ps.ls, func = compute_DESeq2)
compile_3_lvl(test_DESeq2, func = compile_DESeq2) %>% 
  write_tsv(paste0(out_path,'/DESeq2.tsv'))

# EdgeR
test_edgeR <- compute_3_lvl(ps.ls, func = compute_edgeR)
compile_3_lvl(test_edgeR, func = compile_edgeR) %>% 
  write_tsv(paste0(out_path,'/edgeR.tsv'))

# MaAsLin2 ### RAREFIED !!
capture_Maaslin_stdout <- compute_3_lvl(ps_rare.ls, compute_Maaslin2, out_path = out_path)
compile_Maaslin(res_path = paste0(out_path,'/Maaslin2/*/*/*/significant_results.tsv')) %>%
  write_tsv(paste0(out_path,'/Maaslin2.tsv'))

# ZicoSeq ### RAREFIED
test_ZicoSeq <- compute_3_lvl(ps_rare.ls, func = compute_ZicoSeq)
compile_3_lvl(test_ZicoSeq, func = compile_ZicoSeq) %>% 
  write_tsv(paste0(out_path,'/ZicoSeq.tsv'))

ps.ls$Genus <- NULL
# RadEmu
test_radEmu <- compute_3_lvl(ps.ls, func = compute_radEmu)
compile_3_lvl(test_radEmu, func = compile_radEmu) %>% 
  write_tsv(paste0(out_path,'/radEmu.tsv'))

### Parsing results
rbind(
  read_tsv(paste0(out_path,'/Maaslin2.tsv')),
  read_tsv(paste0(out_path,'/AncomBC2.tsv')),
  read_tsv(paste0(out_path,'/corncob.tsv')),
  read_tsv(paste0(out_path,'/edgeR.tsv')), # Too many taxa, needs dealing with !
  read_tsv(paste0(out_path,'/DESEq2.tsv')),
  read_tsv(paste0(out_path,'/radEmu.tsv')),
  read_tsv(paste0(out_path,'/Aldex2.tsv')),
  read_tsv(paste0(out_path,'/ZicoSeq.tsv'))
) %>% dplyr::filter(adj.p < 0.05) %>% 
  write_tsv(paste0(out_path,'/Compiled_DAA_05.tsv'))
