# Are the sourmash benchmark profiles the same with 4.9.4 + branchwater 0.9.15 as with 4.8.11? Part 2 of 2, after
# scripts/benchmark/sourmash_version_check.sh. Three levels, per sample x database:
#   sketch   same md5 of the sample sketch in both versions?
#   genomes  same gather matches (md5), same order, same f_unique_weighted?
#   species  profile from the production tax step + parse_merged_sourmash(): richness, Jaccard, Bray-Curtis,
#            Shannon, fraction classified, top 10 species
# Plus an anchor: the 4.8.11 benchmark summary against the production merged table (data/<ds>/SM_<db>/), which the
# paper uses, for the datasets that have one.
# Output: Out/benchmark/sourmash_version_check.tsv (one row per sample x database) and a verdict per database.
# Run from the repo root: module load StdEnv/2023 r/4.4.0; Rscript scripts/benchmark/compare_sourmash_versions.R
suppressPackageStartupMessages(library(tidyverse))
source('scripts/utils/parse_merged.R')

B <- 'Out/benchmark'
V <- file.path(B, 'version_check')
dbs <- c(RefSeq = 'RefSeq_20250528', GTDB = 'gtdb-rs220-reps')
samples <- read_tsv(file.path(B, 'samples.tsv'), col_names = FALSE, col_types = cols(.default = 'c')) %>%
  transmute(dataset = X1, sample = X2)
tol <- 1e-9 # float rounding

# Sketches ########
sig <- map_dfr(c(old = '4.8.11', new = '4.9.4'), \(v) read_csv(file.path(V, paste0('sigs_', v, '.csv')),
  col_types = cols(.default = 'c')) %>% select(sample = name, md5), .id = 'version') %>%
  pivot_wider(names_from = version, values_from = md5)
sketch <- sig %>% transmute(sample, sketch_same = old == new)

# Genomes (raw gather) ########
# 4.8.11 gather names the match columns name/md5, fastgather match_name/match_md5. Missing or header-only = no match.
read_gather <- function(f) {
  empty <- tibble(match_md5 = character(), match_name = character(), rank = integer(), f_uw = double())
  if (!file.exists(f)) return(empty)
  d <- read_csv(f, col_types = cols(.default = 'c'))
  if (!nrow(d)) return(empty)
  d %>% rename(any_of(c(match_md5 = 'md5', match_name = 'name'))) %>%
    transmute(match_md5, match_name, rank = as.integer(gather_result_rank), f_uw = as.numeric(f_unique_weighted))
}
gather_file <- function(dir, s, db) file.path(dir, paste0(s, '_', db, '_gather.csv'))

genomes <- map_dfr(dbs, \(db) {
  new_dir <- file.path(B, 'sm49', paste0('SM_', db))
  done <- samples$sample[file.exists(gather_file(new_dir, samples$sample, db))] # 4.9.4 gather finished
  map_dfr(done, \(s) {
    o <- read_gather(gather_file(file.path(B, paste0('SM_', db)), s, db))
    n <- read_gather(gather_file(new_dir, s, db))
    j <- full_join(o, n, by = 'match_md5', suffix = c('_old', '_new'))
    tibble(sample = s, genomes_old = nrow(o), genomes_new = nrow(n), genomes_shared = sum(j$match_md5 %in% intersect(o$match_md5, n$match_md5)),
           same_order = nrow(o) == nrow(n) && identical(o$match_md5, n$match_md5),
           max_d_f_unique_weighted = if (nrow(j)) max(abs(coalesce(j$f_uw_old, 0) - coalesce(j$f_uw_new, 0))) else 0)
  })
}, .id = 'database')

# Species profiles ########
summary_file <- function(db, v) file.path(V, paste0(db, '_', v, '.summary.csv.gz'))
# Long table: sample, lineage (7 ranks), count; as in production (parse_merged_sourmash -> species_glom)
profile_long <- function(file) {
  parse_merged_sourmash(file) %>%
    unite('lineage', Kingdom:Species, sep = ';') %>%
    pivot_longer(-lineage, names_to = 'sample', values_to = 'count') %>%
    filter(count > 0)
}
classified <- function(file) { # fraction of the weighted hashes with a species
  read_csv(file, col_types = cols(.default = 'c', f_weighted_at_rank = 'd')) %>%
    filter(rank == 'species', lineage != 'unclassified') %>%
    summarise(f_classified = sum(f_weighted_at_rank), .by = query_name) %>% rename(sample = query_name)
}
shannon <- function(p) -sum(p[p > 0] * log(p[p > 0]))

species <- map_dfr(dbs, \(db) {
  p <- map(c(old = 'old', new = 'new'), \(v) profile_long(summary_file(db, v)) %>%
    mutate(relab = count / sum(count), .by = sample))
  fc <- full_join(classified(summary_file(db, 'old')), classified(summary_file(db, 'new')),
                  by = 'sample', suffix = c('_old', '_new'))
  top10 <- function(d) d %>% slice_max(relab, n = 10, with_ties = FALSE, by = sample) %>%
    summarise(top = list(sort(lineage)), .by = sample)
  j <- full_join(p$old, p$new, by = c('sample', 'lineage'), suffix = c('_old', '_new')) %>%
    mutate(across(c(count_old, count_new, relab_old, relab_new), \(x) coalesce(x, 0)))
  j %>% summarise(
    species_old = sum(count_old > 0), species_new = sum(count_new > 0),
    jaccard = sum(count_old > 0 & count_new > 0) / n(),
    bray_curtis = 0.5 * sum(abs(relab_old - relab_new)),
    max_d_count = max(abs(count_old - count_new)),
    shannon_old = shannon(relab_old), shannon_new = shannon(relab_new),
    only_old_n = sum(count_old > 0 & count_new == 0), only_old_relab = sum(relab_old[count_new == 0]),
    only_new_n = sum(count_new > 0 & count_old == 0), only_new_relab = sum(relab_new[count_old == 0]),
    .by = sample) %>%
    left_join(full_join(top10(p$old), top10(p$new), by = 'sample', suffix = c('_old', '_new')) %>%
      transmute(sample, top10_same = map2_lgl(top_old, top_new, identical)), by = 'sample') %>%
    left_join(fc, by = 'sample')
}, .id = 'database')

# Anchor: 4.8.11 benchmark summary = production merged table? ########
anchor <- map_dfr(dbs, \(db) {
  bench <- read_csv(summary_file(db, 'old'), col_types = cols(.default = 'c', f_weighted_at_rank = 'd')) %>%
    filter(rank == 'species')
  map_dfr(unique(samples$dataset), \(ds) {
    f <- file.path('data', ds, paste0('SM_', db), 'merged_tax_summary.csv.gz')
    if (!file.exists(f)) return(NULL)
    s <- samples$sample[samples$dataset == ds]
    prod <- read_csv(f, col_types = cols(.default = 'c', f_weighted_at_rank = 'd')) %>%
      filter(rank == 'species', query_name %in% s)
    full_join(filter(bench, query_name %in% s), prod, by = c('query_name', 'lineage'), suffix = c('_b', '_p')) %>%
      summarise(anchor_same = all(abs(coalesce(f_weighted_at_rank_b, -1) - coalesce(f_weighted_at_rank_p, -1)) < tol),
                .by = query_name) %>%
      rename(sample = query_name)
  })
}, .id = 'database')

# Table and verdict ########
res <- genomes %>%
  left_join(sketch, by = 'sample') %>%
  left_join(species, by = c('database', 'sample')) %>%
  left_join(anchor, by = c('database', 'sample')) %>%
  mutate(across(c(species_old, species_new, only_old_n, only_new_n), \(x) coalesce(x, 0L)),
         across(c(bray_curtis, only_old_relab, only_new_relab, max_d_count), \(x) coalesce(x, 0)),
         top10_same = coalesce(top10_same, species_old == 0 & species_new == 0)) %>%
  left_join(samples, by = 'sample') %>%
  relocate(database, dataset, sample)
write_tsv(res, file.path(B, 'sourmash_version_check.tsv'))

for (d in names(dbs)) {
  r <- filter(res, database == d)
  identical_ <- with(r, same_order & max_d_f_unique_weighted < tol & bray_curtis < tol)
  equivalent <- with(r, bray_curtis <= 0.05 & top10_same)
  verdict <- if (all(identical_)) 'identical' else if (all(equivalent)) 'equivalent' else 'DIFFERENT'
  cat(sprintf('\n%s: %s (%d samples compared, %d not yet run with 4.9.4)\n', d, verdict, nrow(r), nrow(samples) - nrow(r)))
  cat(sprintf('  sketches identical: %d/%d; genome lists identical: %d; max |d f_unique_weighted|: %.2g\n',
              sum(r$sketch_same, na.rm = TRUE), sum(!is.na(r$sketch_same)), sum(r$same_order), max(r$max_d_f_unique_weighted)))
  cat(sprintf('  Bray-Curtis max %.3g, median %.3g; Jaccard min %.3f; top 10 same in %d; max |d count| %g\n',
              max(r$bray_curtis), median(r$bray_curtis), min(r$jaccard, na.rm = TRUE), sum(r$top10_same), max(r$max_d_count)))
  cat(sprintf('  anchor (4.8.11 benchmark = production table): %d/%d samples\n',
              sum(r$anchor_same, na.rm = TRUE), sum(!is.na(r$anchor_same))))
  if (!all(identical_)) print(r[!identical_, ] %>% select(dataset, sample, genomes_old, genomes_new, genomes_shared,
                                                          bray_curtis, jaccard, only_old_relab, only_new_relab), n = Inf)
}
