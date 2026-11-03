# Readers for the merged tables written by scripts/src/merge_tables.sh, and the per-dataset parsing (meta_parsing()).
# They replace community_functions.R (jorondo1/misc_scripts) parse_MPA() / parse_SM() + lineage joins, which read the
# per-sample files, and return the same tables: one row per taxon, ranks Kingdom..Species (character) and one double
# column per sample, ready for mgx.tools::assemble_phyloseq().
# Needs tidyverse, magrittr, phyloseq, mgx.tools.

# Rank columns from an MPA-style lineage string (k__X|p__X|...|s__X), as parse_MPA() does
mpa_ranks <- function(df, col = 'Taxonomy') {
  tax <- df[[col]]
  df %>%
    dplyr::mutate(
      Kingdom = str_extract(tax, "k__[^|]+") %>% str_remove("k__"),
      Phylum = str_extract(tax, "p__[^|]+") %>% str_remove("p__"),
      Class = str_extract(tax, "c__[^|]+") %>% str_remove("c__"),
      Order = str_extract(tax, "o__[^|]+") %>% str_remove("o__"),
      Family = str_extract(tax, "f__[^|]+") %>% str_remove("f__"),
      Genus = str_extract(tax, "g__[^|]+") %>% str_remove("g__"),
      Species = str_extract(tax, "s__[^|]+") %>% str_remove("s__")) %>%
    dplyr::select(-all_of(col))
}

# Sum genomes (or strains) with the same lineage, as community_functions.R species_glom()
species_glom <- function(abundTable) {
  abundTable %>%
    group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    dplyr::summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)), .groups = "drop")
}

# Merged MPA-style table (merge_tables.sh merged_profile.tsv) -> wide table like parse_MPA()
# tool: 'metaphlan' (merge_metaphlan_tables.py), 'motus3' (merge_tables.sh outer join), 'motus4' (motus merge),
#       'kraken' (combine_mpa.py, header rewritten with sample names)
parse_merged_profile <- function(file, tool = c('metaphlan', 'motus3', 'motus4', 'kraken')) {
  tool <- match.arg(tool)
  lines <- readLines(file)
  header_pattern <- switch(tool,
    metaphlan = '^clade_name\t', motus3 = '^mOTU\t', motus4 = '^mOTU\t', kraken = '^#Classification\t')
  h <- grep(header_pattern, lines)[1]
  if (is.na(h)) stop('no header line in ', file)
  header <- strsplit(lines[h], '\t', fixed = TRUE)[[1]]
  body <- lines[-seq_len(h)] %>% grep('^#', ., invert = TRUE, value = TRUE)
  n_lead <- switch(tool, metaphlan = 1, motus3 = 3, motus4 = 2, kraken = 1)
  lead_names <- switch(tool,
    metaphlan = 'Taxonomy', motus3 = c('mOTU', 'Taxonomy', 'NCBI'), motus4 = c('mOTU', 'Taxonomy'), kraken = 'Taxonomy')

  # read_filename() dropped any per-sample line containing 'Taxonomy' or 'unassigned' (header lines, unassigned mOTUs).
  # The merged rows hold the same text plus the values, so the same patterns drop the same rows.
  body <- body %>%
    grep("Taxonomy", ., invert = TRUE, value = TRUE) %>%
    grep("unassigned", ., invert = TRUE, value = TRUE)

  wide <- read.table(text = body, sep = '\t', header = FALSE, quote = "", comment.char = "",
                     colClasses = 'character', col.names = c(lead_names, header[-seq_len(n_lead)]),
                     check.names = FALSE)

  long <- wide %>%
    tibble() %>%
    pivot_longer(-all_of(lead_names), names_to = 'sample', values_to = 'Abundance') %>%
    dplyr::mutate(Abundance = as.numeric(Abundance)) %>%
    dplyr::filter(Abundance != 0) # absent from that sample's profile (the merge utilities fill with 0)

  # From here on, the steps of parse_MPA()
  long %>%
    dplyr::filter(str_detect(Taxonomy, "s__") &
                    !str_detect(Taxonomy, "t__") &
                    !str_detect(Taxonomy, "UNCLASSIFIED")) %>%
    { if (tool %in% c('motus3', 'motus4')) mutate(., Taxonomy = paste(Taxonomy, mOTU)) else . } %>%
    { # MetaPhlAn: each sample's species sum to 100
      if (tool == 'metaphlan') {
        group_by(., sample) %>% mutate(Abundance = 100 * Abundance / sum(Abundance)) %>% ungroup()
      } else . } %>%
    dplyr::select(sample, Taxonomy, Abundance) %>%
    group_by(Taxonomy, sample) %>%
    dplyr::summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = sample, values_from = Abundance, values_fill = 0) %>%
    mpa_ranks()
}

# sourmash tax metagenome csv_summary of all samples (merge_tables.sh merged_tax_summary.csv.gz) -> wide table like
# parse_SM() + lineage join + species_glom(): species-rank counts = f_weighted_at_rank x total_weighted_hashes
# (= the summed unique_intersect_bp / scaled x average_abund of the species' genomes, rounded as parse_SM() does).
# Matches without lineage (RefSeq plasmid bundles) are in 'unclassified', which is dropped, as parse_SM() drops plasmids.
parse_merged_sourmash <- function(file) {
  read_csv(file, show_col_types = FALSE,
           col_types = cols(.default = 'c', f_weighted_at_rank = 'd', total_weighted_hashes = 'd')) %>%
    dplyr::filter(rank == 'species', lineage != 'unclassified') %>%
    dplyr::transmute(
      sample = query_name,
      lineage,
      count = round(f_weighted_at_rank * total_weighted_hashes, digits = 0)) %>%
    tidyr::separate_wider_delim(lineage, delim = ';', too_few = 'align_start',
      names = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')) %>%
    # parse_GTDB_lineages()/parse_genbank_lineages() strip the rank prefixes (d__, p__, ...)
    dplyr::mutate(across(Kingdom:Species, ~str_remove(., "^[A-Za-z]_+"))) %>%
    tidyr::pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
    species_glom()
}

# sylph: counts.tsv (genome x sample read counts, already rate-corrected by sylph_merge.sh) + taxonomy.tsv (same rows)
parse_sylph <- function(sylph_dir) {
  left_join(
    read_tsv(file.path(sylph_dir, 'counts.tsv'), show_col_types = FALSE),
    read_tsv(file.path(sylph_dir, 'taxonomy.tsv'), col_types = cols(.default = 'c')),
    by = 'genome'
  ) %>%
    dplyr::rename(Kingdom = Domain) %>%
    # genomes without species (no lineage, e.g. the RefSeq plasmid.* bundles) can't be assigned
    dplyr::filter(!is.na(Species)) %>%
    dplyr::select(-genome) %>%
    dplyr::mutate(across(where(is.numeric), as.double))
}

# Database folders of a dataset, by tool. sourmash: the full GTDB rs220 gathers (SM_gtdb_rs220, P19 only) aren't used.
tool_dirs <- function(dsName) {
  dirs <- list.dirs(file.path('data', dsName), recursive = FALSE) %>% basename()
  list(
    motus3 = intersect('MOTUS', dirs),
    motus4 = intersect('MOTUS4', dirs),
    sourmash = grep('^SM_', dirs, value = TRUE) %>% setdiff('SM_gtdb_rs220'),
    sylph = grep('^SYLPH_', dirs, value = TRUE),
    kraken = grep('^KB', dirs, value = TRUE),
    metaphlan = grep('^MPA_', dirs, value = TRUE)
  )
}

# One phyloseq object per database of a dataset, from the merged tables: List$<database>.ps
# Object names as before: MOTUS3, MOTUS4, SM_<db>, SYLPH_<db>, KB<x>, MPA_<db>
meta_parsing <- function(dsName, samData) {
  ps <- list()
  dirs <- tool_dirs(dsName)
  merged <- function(db, f = 'merged_profile.tsv') {
    path <- file.path('data', dsName, db, f)
    if (!file.exists(path)) stop(path, ' not found: run scripts/src/merge_tables.sh ', dsName, ' on the HPC')
    path
  }
  assemble <- function(tbl, ...) assemble_phyloseq(tbl, samData, seqdepth_plot = FALSE, onlySpecies = TRUE, ...)

  for (db in dirs$motus3) {
    message(paste('Parsing', db, '...'))
    ps[['MOTUS3']] <- parse_merged_profile(merged(db), 'motus3') %>% assemble()
  }
  for (db in dirs$motus4) {
    message(paste('Parsing', db, '...'))
    ps[['MOTUS4']] <- parse_merged_profile(merged(db), 'motus4') %>% assemble()
  }
  for (db in dirs$sourmash) {
    message(paste('Parsing', db, '...'))
    ps[[db]] <- parse_merged_sourmash(merged(db, 'merged_tax_summary.csv.gz')) %>% assemble()
  }
  # sylph: dirs without counts.tsv aren't merged yet: skipped
  for (db in dirs$sylph[file.exists(file.path('data', dsName, dirs$sylph, 'counts.tsv'))]) {
    message(paste('Parsing', db, '...'))
    ps[[db]] <- parse_sylph(file.path('data', dsName, db)) %>% species_glom() %>% assemble()
  }
  for (db in dirs$kraken) {
    message(paste('Parsing', db, '...'))
    ps[[db]] <- parse_merged_profile(merged(db), 'kraken') %>% assemble()
  }
  for (db in dirs$metaphlan) {
    message(paste('Parsing', db, '...'))
    ps[[db]] <- parse_merged_profile(merged(db), 'metaphlan') %>%
      assemble(min_sample_size = 0, min_taxa_count = 0) # MetaPhlAn is proportions
  }
  ps
}
