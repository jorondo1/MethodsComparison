# Readers for the merged tables written by scripts/src/merge_tables.sh, and the per-dataset parsing (meta_parsing()).
# They replace community_functions.R (jorondo1/misc_scripts) parse_MPA() / parse_SM() + lineage joins, which read the
# per-sample files, and return the same tables: one row per taxon, ranks Kingdom..Species (character) and one double
# column per sample, ready for mgx.tools::assemble_phyloseq().
# Needs tidyverse, magrittr, phyloseq, mgx.tools.

# Rank columns from an MPA-style lineage string (k__X|p__X|...|s__X), as parse_MPA() does.
# Kingdom also from d__ (GTDB domain, mOTUs 4 and sylph); ranks also separated by ';' (mOTUs 4).
mpa_ranks <- function(df, col = 'Taxonomy') {
  tax <- df[[col]]
  rank <- function(prefix) str_extract(tax, paste0("(?<=(^|[|;])", prefix, "__)[^|;]+"))
  df %>%
    dplyr::mutate(
      Kingdom = coalesce(rank('k'), rank('d')),
      Phylum = rank('p'), Class = rank('c'), Order = rank('o'), Family = rank('f'), Genus = rank('g'),
      Species = str_extract(tax, "(?<=s__)[^|]+")) %>% # species last: may hold ';'-free text with spaces
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

# sylph (sylph_merge.sh), species rows, clade_name = d__X|p__X|...|s__X:
#   abundance = 'sequence'  : merged_species.tsv, estimated read counts (rate-corrected)
#   abundance = 'taxonomic' : merged_species_tax.tsv, taxonomic abundance (% of genomes; rescaled in prokaryotes_only())
# Genomes without lineage (NO_TAXONOMY, e.g. the RefSeq plasmid.* bundles) have no species row, so they're left out.
parse_sylph <- function(sylph_dir, abundance = c('sequence', 'taxonomic')) {
  abundance <- match.arg(abundance)
  f <- c(sequence = 'merged_species.tsv', taxonomic = 'merged_species_tax.tsv')[[abundance]]
  read_tsv(file.path(sylph_dir, f), show_col_types = FALSE,
           col_types = cols(clade_name = 'c', .default = 'd')) %>%
    dplyr::mutate(clade_name = str_replace(clade_name, '^d__', 'k__')) %>%
    mpa_ranks('clade_name') %>%
    dplyr::relocate(Kingdom:Species)
}

# Bacteria and Archaea only (the study compares prokaryotic profiles; RefSeq and GTDB databases made comparable).
# Drops Eukaryota (human, fungi, Blastocystis...), viruses (Viruses, or no domain in RefSeq/NCBI lineages: phages, phiX)
# and taxa without a domain. Plasmid bundles and NO_TAXONOMY genomes have no species and are already gone.
# Kraken GTDB lineages have no domain level (they start at p__): GTDB holds Bacteria and Archaea only (all 159 phyla
# found in the KB*_GTDB tables are GTDB r220 bacterial or archaeal phyla), so they're kept whole.
# rescale: proportions (MetaPhlAn, sylph taxonomic abundance) are rescaled to 100 per sample after the filter.
prokaryotes_only <- function(tbl, db, rescale = FALSE) {
  if (str_detect(db, '_GTDB$')) {
    out <- tbl
  } else {
    out <- dplyr::filter(tbl, Kingdom %in% c('Bacteria', 'Archaea'))
  }
  if (rescale) out <- dplyr::mutate(out, across(where(is.numeric), \(x) if (sum(x) > 0) 100 * x / sum(x) else x))
  out
}

# Database folders of a dataset, by tool. Not used: the full GTDB rs220 gathers (SM_gtdb_rs220, P19 only) and Maple's
# sylph fungi RefSeq and GlobDB r232 profiles (side tests).
tool_dirs <- function(dsName) {
  dirs <- list.dirs(file.path('data', dsName), recursive = FALSE) %>% basename()
  list(
    motus3 = intersect('MOTUS', dirs),
    motus4 = intersect('MOTUS4', dirs),
    sourmash = grep('^SM_', dirs, value = TRUE) %>% setdiff('SM_gtdb_rs220'),
    sylph = grep('^SYLPH_', dirs, value = TRUE) %>% setdiff(c('SYLPH_fungi_refseq', 'SYLPH_glob_r232')), # Maple side tests
    kraken = grep('^KB', dirs, value = TRUE),
    metaphlan = grep('^MPA_', dirs, value = TRUE)
  )
}

# One phyloseq object per database of a dataset, from the merged tables: List$<database>.ps
# Object names: MOTUS3, MOTUS4, SM_<db>, SYLPH_<db>_Seq, SYLPH_<db>_Tax, KB<x>, MPA_<db>
meta_parsing <- function(dsName, samData) {
  ps <- list()
  dirs <- tool_dirs(dsName)
  merged <- function(db, f = 'merged_profile.tsv') {
    path <- file.path('data', dsName, db, f)
    if (!file.exists(path)) stop(path, ' not found: run scripts/src/merge_tables.sh (sylph: sylph_merge.sh) ', dsName, ' on the HPC')
    path
  }
  assemble <- function(tbl, ...) assemble_phyloseq(tbl, samData, seqdepth_plot = FALSE, onlySpecies = TRUE, ...)

  for (db in dirs$motus3) {
    message(paste('Parsing', db, '...'))
    ps[['MOTUS3']] <- parse_merged_profile(merged(db), 'motus3') %>% prokaryotes_only('MOTUS3') %>% assemble()
  }
  for (db in dirs$motus4) {
    message(paste('Parsing', db, '...'))
    ps[['MOTUS4']] <- parse_merged_profile(merged(db), 'motus4') %>% prokaryotes_only('MOTUS4') %>% assemble()
  }
  for (db in dirs$sourmash) {
    message(paste('Parsing', db, '...'))
    ps[[db]] <- parse_merged_sourmash(merged(db, 'merged_tax_summary.csv.gz')) %>% prokaryotes_only(db) %>% assemble()
  }
  # sylph: two objects per database, <db>_Seq (read counts) and <db>_Tax (taxonomic abundance, not rarefied).
  # Dirs without merged_species.tsv aren't merged yet: skipped
  for (db in dirs$sylph[file.exists(file.path('data', dsName, dirs$sylph, 'merged_species.tsv'))]) {
    message(paste('Parsing', db, '...'))
    ps[[paste0(db, '_Seq')]] <- parse_sylph(file.path('data', dsName, db), 'sequence') %>% species_glom() %>%
      prokaryotes_only(db) %>% assemble()
    merged(db, 'merged_species_tax.tsv') # stops if missing: sylph_merge.sh --tax-only <dataset> <db name>
    ps[[paste0(db, '_Tax')]] <- parse_sylph(file.path('data', dsName, db), 'taxonomic') %>%
      species_glom() %>% prokaryotes_only(db, rescale = TRUE) %>%
      assemble(min_sample_size = 0, min_taxa_count = 0) # proportions
  }
  for (db in dirs$kraken) {
    message(paste('Parsing', db, '...'))
    ps[[db]] <- parse_merged_profile(merged(db), 'kraken') %>% prokaryotes_only(db) %>% assemble()
  }
  for (db in dirs$metaphlan) {
    message(paste('Parsing', db, '...'))
    ps[[db]] <- parse_merged_profile(merged(db), 'metaphlan') %>% prokaryotes_only(db, rescale = TRUE) %>%
      assemble(min_sample_size = 0, min_taxa_count = 0) # MetaPhlAn is proportions
  }
  ps
}

# Dataset-specific changes to meta_parsing()'s list (formerly inline in 0_proc_data.R).
# sam_data: load_sample_data()'s list (Moss's MAGs object is in it).
dataset_fixes <- function(dsName, ps, sam_data) {
  if (dsName == 'Moss') {
    ps[['SM_gtdb-rs214-rep_MAGs']] <- sam_data$moss.ps # sourmash GTDB 214 + the Moss MAGs, from the external RDS
    ps$MPA_db2022 <- NULL
    ps$MPA_db2023 <- NULL
    ps$KB90 <- NULL
    ps$MOTUS3 <- NULL # NULL anyway: no Moss sample reaches min_sample_size
  }
  if (dsName == 'Maple') {
    # Maple's GTDB rs220 representatives gather dir is spelled SM_gtdb-rs220-reps: same database as the other datasets' SM_gtdb-rs220-rep
    names(ps) <- str_replace(names(ps), '^SM_gtdb-rs220-reps$', 'SM_gtdb-rs220-rep')
  }
  ps
}
