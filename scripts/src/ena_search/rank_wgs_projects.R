# Rank ENA metagenomic WGS studies as candidate public datasets:
# non-human, >= min_samples samples, and likely two-group design.
# Input: Out/ena_candidates/runs_2023.tsv from find_wgs_projects.sh
# Run from the project root: Rscript scripts/src/ena_search/rank_wgs_projects.R

library(pacman)
p_load(tidyverse, xml2, httr2)

out_dir <- 'Out/ena_candidates'
runs_file <- file.path(out_dir, 'runs_2023.tsv')
min_samples <- 50
n_shortlist <- 150        # studies whose sample attributes are fetched (stage 3)
max_attr_samples <- 300   # samples per study fetched for attribute checks
xml_batch_size <- 100

set.seed(1)

#### Stage 2: per-study summary ####

runs <- read_tsv(runs_file,
                 col_types = cols(.default = 'c', base_count = 'd'),
                 quote = '', na = '')

mode_chr <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  names(which.max(table(x)))
}

samples <- runs %>%
  group_by(study_accession, sample_accession) %>%
  summarise(scientific_name = first(scientific_name),
            sample_alias = first(sample_alias),
            sample_title = first(sample_title),
            bases = sum(base_count, na.rm = TRUE),
            .groups = 'drop')

studies <- runs %>%
  group_by(study_accession) %>%
  summarise(study_title = first(study_title),
            first_public = min(first_public),
            n_runs = n(),
            platform = mode_chr(instrument_platform),
            layout = mode_chr(library_layout),
            .groups = 'drop') %>%
  left_join(samples %>%
              group_by(study_accession) %>%
              summarise(n_samples = n(),
                        scientific_name = mode_chr(scientific_name),
                        n_sci_names = n_distinct(scientific_name),
                        median_gbp = round(median(bases) / 1e9, 2),
                        .groups = 'drop'),
            by = 'study_accession') %>%
  filter(!str_detect(coalesce(study_title, ''),
                     regex('amplicon|16S|Should be deleted', ignore_case = TRUE)))

#### Habitat category ####

# Metagenome names that don't say where the sample came from; use the title instead
generic_names <- regex(paste0(
  '^(metagenome|(gut|feces|fecal|oral|skin|vaginal|blood|lung|respiratory tract|',
  'upper respiratory tract|milk|urine|synthetic|viral|mixed culture) metagenome)$'),
  ignore_case = TRUE)

human_title <- paste0(
  '\\bhumans?\\b|homo sapiens|patient|\\bchild|infant|neonat|cohort|\\bwomen\\b|\\bmen\\b|',
  'volunteer|participant|clinical|elderly|pregnan|\\bpeople\\b|individuals')

# Ordered: first match wins
env_patterns <- c(
  Extreme = 'hot spring|hypersaline|salt lake|saline lake|acid mine|hydrothermal|halophil|thermophil|hypolith|cryoconite',
  Sediment = 'sediment|\\bmud\\b',
  Marine = 'marine|seawater|sea water|ocean|estuar|coral|sponge|seagrass|coastal|\\bseas?\\b|\\bbay\\b',
  Freshwater = 'freshwater|fresh water|lake|river|aquatic|\\bponds?\\b|groundwater|ground water|drinking water|\\bstreams?\\b|reservoir|wetland',
  Engineered = 'wastewater|sewage|sludge|bioreactor|digester|biofilm|indoor|surface metagenome|urban|\\bair\\b|dust|clothing|aquaculture|manure|built environment|treatment plant',
  Food = 'food|ferment|dairy|cheese|kimchi|wine|beer|sourdough|\\bmeat\\b|vegetable',
  Animal = paste0(
    'mouse|mice|murine|\\brats?\\b|rodent|porcine|\\bpigs?\\b|piglet|swine|bovine|cattle|\\bcows?\\b|',
    'calf|calves|rumen|ruminant|chicken|poultry|broiler|\\bhens?\\b|duck|goose|fish|zebrafish|',
    'salmon|insect|\\bbees?\\b|honeybee|bumblebee|\\bants?\\b|termite|mosquito|beetle|larva|',
    'mealworm|silkworm|sheep|\\bgoats?\\b|canine|\\bdogs?\\b|feline|\\bcats?\\b|horse|equine|',
    'primate|monkey|macaque|invertebrate|mollus|oyster|mussel|shrimp|crab|\\bbirds?\\b|avian|',
    '\\bbats?\\b|deer|\\byaks?\\b|camel|rabbit|wildlife|animal|\\bticks?\\b|nematode|earthworm|',
    'amphibian|\\bfrogs?\\b|reptile|turtle|snake|lizard|whale|dolphin|\\bseals?\\b|panda|vertebrate'),
  Plant = 'plant|rhizosphere|\\broots?\\b|\\bleaf|leaves|phyllosphere|endophyt|\\bseeds?\\b|\\brice\\b|wheat|maize|\\bcorn\\b|\\bmoss|\\btrees?\\b|\\bcrops?\\b|lichen',
  Soil = 'soil|permafrost|\\bpeat|compost|desert|rhizosheath'
)

classify <- function(x, patterns) {
  out <- rep(NA_character_, length(x))
  for (h in names(patterns)) {
    hit <- is.na(out) & coalesce(str_detect(x, regex(patterns[[h]], ignore_case = TRUE)), FALSE)
    out[hit] <- h
  }
  out
}

studies <- studies %>%
  mutate(
    habitat_sci = case_when(
      str_detect(scientific_name, regex('^human', ignore_case = TRUE)) ~ 'Human',
      coalesce(str_detect(scientific_name, generic_names), TRUE) ~ NA_character_,
      TRUE ~ classify(scientific_name, env_patterns)),
    habitat = coalesce(habitat_sci,
                       classify(study_title, c(Human = human_title, env_patterns)),
                       'Uncertain')) %>%
  select(-habitat_sci)

#### Two-group split detection from sample names ####

no_split <- tibble(balance = 0, score = 0, desc = NA_character_)

# Two most frequent values cover >= cover of n samples, and the smaller one >= min_share
value_split <- function(x, n = length(x), cover = 0.9, min_share = 0.2) {
  x <- x[!is.na(x) & x != '']
  if (!length(x)) return(no_split)
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) < 2) return(no_split)
  top <- tab[1:2]
  if (sum(top) / n < cover || top[2] / n < min_share) return(no_split)
  bal <- as.numeric(top[2] / top[1])
  tibble(balance = bal, score = bal,
         desc = sprintf('%s(%d) | %s(%d)', names(top)[1], top[1], names(top)[2], top[2]))
}

# Groups encoded as words in sample names, e.g. rhizo_12 vs bulk_12, or Abx present/absent.
# Single-token (present vs absent) splits are noisier, so their score is down-weighted.
token_split <- function(labels, max_tokens = 30) {
  n <- length(labels)
  toks <- map(str_extract_all(str_to_lower(coalesce(labels, '')), '[a-z]{2,}'), unique)
  share <- table(unlist(toks)) / n
  cand <- names(share)[share >= 0.15 & share <= 0.85]
  if (!length(cand)) return(no_split)
  cand <- head(cand[order(-share[cand])], max_tokens)
  M <- vapply(cand, function(t) map_lgl(toks, ~ t %in% .x), logical(n))
  counts <- colSums(M)
  best <- no_split

  if (length(cand) >= 2) {
    co <- crossprod(M)
    pairs <- which(co == 0 & upper.tri(co), arr.ind = TRUE)
    if (nrow(pairs)) {
      a <- counts[pairs[, 1]]
      b <- counts[pairs[, 2]]
      bal <- ifelse((a + b) / n >= 0.8, pmin(a, b) / pmax(a, b), 0)
      i <- which.max(bal)
      if (bal[i] > 0) {
        best <- tibble(balance = bal[i], score = bal[i],
                       desc = sprintf('%s(%d) | %s(%d)',
                                      cand[pairs[i, 1]], a[i], cand[pairs[i, 2]], b[i]))
      }
    }
  }

  single <- counts / n >= 0.2 & counts / n <= 0.8
  if (any(single)) {
    bal <- ifelse(single, pmin(counts, n - counts) / pmax(counts, n - counts), 0)
    i <- which.max(bal)
    if (0.75 * bal[i] > best$score) {
      best <- tibble(balance = bal[i], score = 0.75 * bal[i],
                     desc = sprintf('%s(%d) | other(%d)', cand[i], counts[i], n - counts[i]))
    }
  }
  best
}

name_split <- function(df) {
  bind_rows(
    alias = token_split(df$sample_alias),
    title = token_split(df$sample_title),
    sci_name = value_split(df$scientific_name, cover = 0.9, min_share = 0.15),
    .id = 'source') %>%
    slice_max(score, n = 1, with_ties = FALSE) %>%
    transmute(name_score = score,
              name_split = if_else(score > 0, paste0(source, ': ', desc), NA_character_))
}

candidates <- studies %>%
  filter(habitat != 'Human', n_samples >= min_samples)

name_splits <- samples %>%
  semi_join(candidates, by = 'study_accession') %>%
  group_by(study_accession) %>%
  group_modify(~ name_split(.x)) %>%
  ungroup()

title_regex <- regex(paste0(
  'versus|\\bvs\\b|compar|between|across|differen|gradient|contrast|effects? of|',
  'response|impact|season|depth|habitat|species'), ignore_case = TRUE)

candidates <- candidates %>%
  left_join(name_splits, by = 'study_accession') %>%
  mutate(title_hits = map_chr(str_extract_all(str_to_lower(coalesce(study_title, '')), title_regex),
                              ~ paste(unique(.x), collapse = ',')),
         title_score = pmin(str_count(title_hits, '[^,]+'), 2) / 2,
         name_score = coalesce(name_score, 0),
         prelim_score = 2 * name_score + title_score)

#### Stage 3: sample attributes for the shortlist ####

boilerplate_tags <- regex(paste0(
  '^(ENA-|INSDC|NCBI|BioSampleModel|SRA accession|collection[ _]date|lat[ _]lon|',
  'geographic location \\((latitude|longitude)\\)|organism$|title$|description$|',
  'sample[ _]name$|submitter|ngdc|alias$|External Id|ENA checklist)|_id$|accession'),
  ignore_case = TRUE)

fetch_sample_attrs <- function(accs) {
  url <- paste0('https://www.ebi.ac.uk/ena/browser/api/xml/', paste(accs, collapse = ','))
  tryCatch({
    doc <- request(url) %>% req_retry(max_tries = 3) %>% req_perform() %>% resp_body_xml()
    attrs <- xml_find_all(doc, '//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')
    tibble(sample_accession = xml_attr(xml_find_first(attrs, './ancestor::SAMPLE'), 'accession'),
           tag = xml_text(xml_find_first(attrs, './TAG')),
           value = xml_text(xml_find_first(attrs, './VALUE')))
  }, error = function(e) {
    warning('Batch starting with ', accs[1], ' failed: ', conditionMessage(e))
    tibble(sample_accession = character(), tag = character(), value = character())
  })
}

shortlist_samples <- candidates %>%
  slice_max(prelim_score, n = n_shortlist, with_ties = FALSE) %>%
  select(study_accession) %>%
  inner_join(samples, by = 'study_accession') %>%
  group_by(study_accession) %>%
  slice_sample(n = max_attr_samples) %>%
  mutate(n_fetched = n()) %>%
  ungroup() %>%
  select(study_accession, sample_accession, n_fetched)

acc_batches <- unique(shortlist_samples$sample_accession) %>%
  split(ceiling(seq_along(.) / xml_batch_size))

message('Fetching sample attributes: ', length(acc_batches), ' batches')
attrs <- imap_dfr(acc_batches, function(accs, i) {
  if (as.integer(i) %% 20 == 0) message('  batch ', i, '/', length(acc_batches))
  Sys.sleep(0.2)
  fetch_sample_attrs(accs)
}) %>%
  filter(!str_detect(tag, boilerplate_tags)) %>%
  inner_join(shortlist_samples, by = 'sample_accession', relationship = 'many-to-many')

# Human studies hidden behind generic metagenome names
human_hosts <- attrs %>%
  filter(str_detect(tag, regex('^host([ _]scientific[ _]name)?$', ignore_case = TRUE))) %>%
  group_by(study_accession) %>%
  summarise(human_share = sum(str_detect(value, regex('homo sapiens|^human$', ignore_case = TRUE))) /
              first(n_fetched)) %>%
  filter(human_share >= 0.5) %>%
  pull(study_accession)

attr_splits <- attrs %>%
  group_by(study_accession, tag) %>%
  summarise(split = list(value_split(value, n = first(n_fetched))), .groups = 'drop') %>%
  unnest(split) %>%
  filter(balance > 0) %>%
  group_by(study_accession) %>%
  slice_max(balance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(study_accession, attr_balance = balance, attr_split = paste0(tag, ': ', desc))

#### Final ranking ####

studies <- studies %>%
  mutate(habitat = if_else(study_accession %in% human_hosts, 'Human', habitat))

ranked <- candidates %>%
  filter(!study_accession %in% human_hosts) %>%
  left_join(attr_splits, by = 'study_accession') %>%
  mutate(attr_balance = coalesce(attr_balance, 0),
         design_score = round(3 * attr_balance + 2 * name_score + title_score, 3),
         ena_url = paste0('https://www.ebi.ac.uk/ena/browser/view/', study_accession)) %>%
  arrange(desc(design_score), desc(n_samples)) %>%
  select(study_accession, ena_url, habitat, scientific_name, study_title, n_samples, n_runs,
         median_gbp, platform, layout, design_score, attr_split, name_split, title_hits,
         first_public)

studies %>%
  arrange(desc(n_samples)) %>%
  write_tsv(file.path(out_dir, 'wgs_projects_2023.tsv'), na = '')
write_tsv(ranked, file.path(out_dir, 'wgs_candidates_ranked.tsv'), na = '')

message('Studies: ', nrow(studies), ' | with >= ', min_samples, ' samples: ',
        sum(studies$n_samples >= min_samples), ' | ranked (non-human): ', nrow(ranked))
print(count(studies, habitat, sort = TRUE), n = Inf)
print(count(ranked, habitat, sort = TRUE), n = Inf)
