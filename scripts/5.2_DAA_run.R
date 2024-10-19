library(pacman)
p_load(magrittr, tidyverse, phyloseq,
       furrr, purrr, parallel,
       ANCOMBC, ALDEx2, radEmu, edgeR, DESeq2, Maaslin2, GUniFrac)

source('scripts/myFunctions.R')
source('scripts/5.1_DAA_fun.R')

ps.ls <- read_rds('Out/ps_full.ls.rds')
ps_rare.ls <- read_rds('Out/ps_rare.ls.rds')

out_path <- 'Out/DAA_Moss'
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

ncores <- detectCores() -1

# Aldex2
test_aldex <- compute_3_lvl(ps.ls, func = compute_aldex)
compile_3_lvl(test_aldex, func = compile_aldex) %>% 
  write_tsv(paste0(out_path,'/Aldex2.tsv'))

# ANCOM-BC2
test_ancombc2 <- compute_3_lvl(ps.ls, func = compute_ancombc2)
compile_3_lvl(test_ancombc2, func = compile_ancombc2) %>% 
  write_tsv(paste0(out_path,'/AncomBC2.tsv'))

# RadEmu
test_radEmu <- compute_3_lvl(ps.ls, func = compute_radEmu)
compile_3_lvl(test_radEmu, func = compile_radEmu) %>% 
  write_tsv(paste0(out_path,'/radEmu.tsv'))

# MaAsLin2 ### RAREFIED !!
capture_Maaslin_stdout <- compute_3_lvl(ps_rare.ls, compute_Maaslin2)
compile_Maaslin(res_path = paste0(out_path,'/Maaslin2/*/*/*/significant_results.tsv')) %>%
  write_tsv(paste0(out_path,'/Maaslin2.tsv'))

# EdgeR
test_edgeR <- compute_3_lvl(ps.ls, func = compute_edgeR)
compile_3_lvl(test_edgeR, func = compile_edgeR) %>% 
  write_tsv(paste0(out_path,'/edgeR.tsv'))

# DESeq2 !Note: L2FC transformed to L10FC for scale
test_DESeq2 <- compute_3_lvl(ps.ls, func = compute_DESeq2)
compile_3_lvl(test_DESeq2, func = compile_DESeq2) %>% 
  write_tsv(paste0(out_path,'/DESeq2.tsv'))

# ZicoSeq ### RAREFIED
test_ZicoSeq <- compute_3_lvl(ps_rare.ls, func = compute_ZicoSeq)
compile_3_lvl(test_ZicoSeq, func = compile_ZicoSeq) %>% 
   write_tsv(paste0(out_path,'/ZicoSeq.tsv'))

### Parsing results
rbind(
  read_tsv(paste0(out_path,'/Maaslin2.tsv')),
  read_tsv(paste0(out_path,'/AncomBC2.tsv')),
  read_tsv(paste0(out_path,'/edgeR.tsv')), # Too many taxa, needs dealing with !
  read_tsv(paste0(out_path,'/DESEq2.tsv')),
  read_tsv(paste0(out_path,'/radEmu.tsv')),
  read_tsv(paste0(out_path,'/Aldex2.tsv')),
  read_tsv(paste0(out_path,'/ZicoSeq.tsv'))
) %>% dplyr::filter(adj.p < 0.05) %>% 
  write_tsv(paste0(out_path,'/Compiled_DAA_05.tsv'))
