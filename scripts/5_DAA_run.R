library(pacman)
p_load(magrittr, tidyverse, phyloseq,
       furrr, purrr, parallel,
       ANCOMBC, ALDEx2, radEmu, edgeR, DESeq2, Maaslin2)

source('scripts/myFunctions.R')
source('scripts/5_DAA_fun.R')

ps.ls <- read_rds('Out/ps_full.ls.rds')

out_path <- 'Out/DAA'
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

ncores <- detectCores() -1

# Aldex2
test_aldex <- compute_3_lvl(ps.ls, func = compute_aldex)
compile_3_lvl(test_aldex, func = compile_aldex) %>% 
  write_tsv('Out/DAA/Aldex2.tsv')

# ANCOM-BC2
test_ancombc2 <- compute_3_lvl(ps.ls, func = compute_ancombc2)
compile_3_lvl(test_ancombc2, func = compile_ancombc2) %>%
  write_tsv('Out/DAA/AncomBC2.tsv')

# RadEmu
test_radEmu <- compute_3_lvl(ps.ls, func = compute_radEmu)
compile_3_lvl(test_radEmu, func = compile_radEmu) %>% 
  write_tsv('Out/DAA/radEmu.tsv')

# MaAsLin2
capture_Maaslin_stdout <- compute_3_lvl(ps.ls, compute_Maaslin2)
compile_Maaslin(res_path = 'Out/DAA/Maaslin2/*/*/*/significant_results.tsv') %>% 
  write_tsv('Out/DAA/Maaslin2.tsv')

# EdgeR
test_edgeR <- compute_3_lvl(ps.ls, func = compute_edgeR)
compile_3_lvl(test_edgeR, func = compile_edgeR) %>% 
  write_tsv('Out/DAA/edgeR.tsv')

# DESeq2
# ...