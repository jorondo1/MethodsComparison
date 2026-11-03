# Merged tables: format checks (phase A)

2026-11-03, branch `merged_tables`, HPC node cv3401 (salloc). Test samples: P19_Saliva GQ10, GQ11, GQ12 (plus AD_Skin for trimmed MetaPhlAn profiles). Outputs in `tmp/merge_tests/<tool>/` (gitignored). Reference logic: `misc_scripts/community_functions.R` at `b86aac6`.

## Summary

| Tool | Utility | Works as is? | Values vs per-sample files |
|---|---|---|---|
| MetaPhlAn 4.2.4 | `merge_metaphlan_tables.py` (metaphlan.4.2.4.sif) | yes | identical (max diff 2e-15, float printing) |
| mOTUs 4.0.4 | `motus merge` (mOTUs_v4.0.4.sif) | yes | identical |
| mOTUs 3.1.0 | `motus merge` (motus3_env) | **no**: refuses our profiles | needs our own outer join |
| Kraken2/Bracken | KrakenTools 1.2 `combine_mpa.py` | yes, but sample names must be set | identical |
| sourmash 4.8.11 | `sourmash tax metagenome -F csv_summary` | per sample only, then concatenate | identical for 4 DBs; 2 DBs need a derived lineage file |

## 1. MetaPhlAn (`merge_metaphlan_tables.py`)

- Command: `singularity exec -B $ILAFORES metaphlan.4.2.4.sif merge_metaphlan_tables.py -o out.tsv <profiles>`.
- Output: line 1 is `#<db version>`, line 2 is `clade_name<TAB><sample>…`, then one row per clade. The NCBI column is dropped.
- Sample columns: the file name with `_profile.txt` stripped (GQ10, SRR20002340), the same as `read_filename()`.
- Full profiles (UNCLASSIFIED + all ranks, e.g. P19_Saliva MPA_db2025) and trimmed ones (`s__` lines only, e.g. AD_Skin) merge together. Absent clades are written as `0.0`.
- Species rows (`s__`, no `t__`) match the raw files: 446/446 values, max abs diff 1.8e-15.
- The reader must keep doing: s__ without t__, drop UNCLASSIFIED, renormalise each sample to 100.

## 2. mOTUs

**4.0.4:** `singularity exec -e -B $ILAFORES mOTUs_v4.0.4.sif motus merge -i a b c -o out.tsv`.
- Loads the DB bundled in the container (~15 s), so no `-db` needed.
- Output: `#tool_version=…` line, then `mOTU<TAB>Taxonomy<TAB><samples>`. mOTU IDs and the full `d__…;s__…` taxonomy are kept. Values stay integer counts, absent = `0`. The `mOTUv4.0_unassigned` row is kept (the reader drops it).
- All 629 non-zero values identical to the per-sample files.

**3.1.0:** `motus merge -db $ILAFORES/ref_dbs/motus3/db_mOTU -i a,b,c` fails: `The taxa ids are different`.
- `motus3_SLURM.sh` pipes the profile through `awk '$4 != 0'`, so each profile lists a different set of mOTUs. `motus merge` (v3) requires identical row sets.
- Plan: a small awk outer join on (mOTU, consensus_taxonomy, NCBI_tax_id), writing the same layout as `motus merge`.
- The `#` header lines are format version 3.1.0 in every dataset.

**Empty profiles (mOTUs 3 only):** PD has 110 empty `MOTUS/*_profile.txt` (all dated 2025-04-10), AD_Skin has 1. `read_filename()` returns 0 rows for them, so those samples are absent from today's MOTUS3 objects. The merge will skip them and list them, so the samples stay absent.

## 3. Kraken (`combine_mpa.py`)

- Command: `python3 $ILAFORES/programs/KrakenTools-1.2/combine_mpa.py -i <files> -o out.tsv`.
- Sample names come from a `#` header line in each input (last tab field). Our `*_bracken_S.MPA.TXT` files have none, so columns are `Sample #1…#N` in input order. The merge script will rewrite the header from the input list (sample = file name minus `_bracken_S.MPA.TXT`).
- Absent taxa are written as `0`. KB10 (1977 values) and KB10_GTDB (6328 values) are identical to the raw files.
- No duplicate taxa within a file (300 random files checked). Every line has 2 fields.
- KB*_GTDB lineages start at `p__` (no `k__`), so Kingdom is NA today too. That doesn't reach the objects: `assemble_phyloseq(onlySpecies = TRUE)` keeps only the Species rank.

## 4. sourmash (`sourmash tax metagenome`)

- **`csv_summary` takes one query at a time.** With several gathers: `Can only output multi-query result formats: skipping csv_summary`. Plan: run it once per sample and concatenate into one long table (`query_name, rank, fraction, lineage, f_weighted_at_rank, bp_match_at_rank, total_weighted_hashes`), one file per dataset × DB.
- **Values:** at the species rank, `f_weighted_at_rank × total_weighted_hashes` equals R's summed `round(unique_intersect_bp / scaled × average_abund)` exactly: the weighted hashes are integers, so R's rounding is a no-op.
- **Unclassified:** each rank has an `unclassified` row. At superkingdom, `1 − fraction(unclassified)` is the classified fraction, so it could replace the sourmash classification-rate block.
- **Lineage format:** GTDB lineages keep their `d__`/`p__` prefixes. NCBI ones (genbank, RefSeq) have none. The reader strips `^[A-Za-z]_+` as today.
- **query_name:** plain sample IDs in every SM dir (no `_clean`).

| DB dir | Lineage file for `sourmash tax` | Result vs R (3 samples) |
|---|---|---|
| SM_genbank-2022.03 | the 5 `genbank-2022.03-*.lineages.csv` | identical |
| SM_gtdb-rs214-full / -rep | `gtdb-rs214.lineages.csv` | identical |
| SM_gtdb-rs220-rep | `gtdb-rs220.lineages.reps.csv` | identical |
| SM_gtdb_rs220 (P19 only) | `gtdb-rs220-full.lineages.csv` | ~70 genomes missing (see below) |
| SM_RefSeq_20250528 | `RefSeq_20250528.lineages.csv` | can't be read (see below) |

**SM_gtdb_rs220 (P19_Gut, P19_Saliva):** this is a gather against the *full* rs220 DB (`gtdb-rs220-k31.zip`, older `name` column). It's not a copy of SM_gtdb-rs220-rep.
- The gathers name `GCA_x` where the lineage file lists the `GCF_x` twin. R joins on the bare 9-digit number, so it doesn't notice, but `sourmash tax` reports them missing.
- Fix: a lineage copy that adds the missing GCA/GCF twin of each ident. With it the 3 samples are identical to R. For phase B, the derived file would be built by the merge script into `$ILAFORES/ref_dbs/sourmash_db/derived/`.

**SM_RefSeq_20250528: today's R objects double-count.**
- `RefSeq_20250528.lineages.csv` has no header, and it lists GCA and GCF twins and several versions of the same accession (112 unversioned IDs with conflicting lineages). `sourmash tax` refuses it (`multiple lineages for identifier GCA_000814125`).
- R's join on the bare number matches most genomes 2+ times: `left_join` duplicates their counts. For GQ10/11/12 the gather counts are 1,925,201 / 220,658 / 217,289, and after R's join 4,148,864 / 465,516 / 469,464 (+111–116 %).
- The stale HPC `Out/_Rdata/ps_raw.ls.RDS` (June 2025) holds 4,142,054 / 446,598 / 453,276 for these samples, so it has the same inflation. The Mac's current objects should be checked too.
- In nearly every case the duplicate rows have the same species, so relative abundances barely change. Absolute counts, rarefaction depth and count-based alpha indices are affected. At least 1 genome per sample is split across 2 species.
- Fix tested: build the sourmash lineage file from the sylph-tax table (`$ILAFORES/ref_dbs/sylph_db/RefSeq_20250528_taxonomy.tsv.gz`: genome file name → lineage, one row per genome), with ident = file name minus `.fna.gz`, and run `sourmash tax` with `--keep-full-identifiers --keep-identifier-versions`.
  - Species totals then equal the gather totals exactly (163/112/118 species).
  - The plasmid bundles become `unclassified`; R drops them too.
  - This is the same taxonomy sylph uses.
- **So the merged sourmash RefSeq objects can't equal today's.** They'll equal the gather totals instead. Decision needed.

## 5. R on the HPC

- `module load StdEnv/2023 r/4.4.0` gives R 4.4.0 with `~/R/x86_64-pc-linux-gnu-library/4.4`. phyloseq, tidyverse and mgx.tools 0.1.0 load. The validation can run here.
- With `onlySpecies = TRUE`, `assemble_phyloseq` keeps only the Species rank. The mOTUs 4 `d__…;` taxonomy (which `parse_MPA`'s `k__[^|]+` extraction mangles) therefore has no effect on the objects.
- The HPC `ps_raw.ls.RDS` predates MOTUS4/MPA_db2025. The comparison has to rebuild the reference from the per-sample files (pinned `b86aac6`), as planned.

## Other things noticed

- Maple has `SM_gtdb-rs220-reps` (with an `s`), not `SM_gtdb-rs220-rep`. `0_config.R` has no entry for it, and none for P19's `SM_gtdb_rs220` either.
- `meta_parsing()` parses every `SM_*` dir, so P19's `SM_gtdb_rs220` (full rs220) is in the objects next to `SM_gtdb-rs220-rep`.

# Results: merged tables vs per-sample parsing (phases B–D)

2026-11-03. Merged with `scripts/src/merge_tables.sh -f --tools motus,mpa,kraken,sourmash` (all 9 datasets of `scripts/src/datasets.tsv`; sylph tables already current). Compared with `scripts/compare_merged.R`: reference = the former `meta_parsing()` blocks on the per-sample files, with `community_functions.R` pinned at `b86aac6`. Full table: `Out/merged_tables/comparison.tsv` (HPC).

| Status | Objects | Which |
|---|---|---|
| identical | 121 | every mOTUs 3/4, Kraken, sylph and non-RefSeq sourmash object; 3 MPA objects |
| equal (float) | 24 | MPA: max relative otu difference 9.7e-13, Shannon 6e-15, richness 0 (merge_metaphlan_tables.py prints the values with different digits) |
| both NULL | 3 | Moss MOTUS3, Bee MOTUS3, Bee MOTUS4: no sample reaches `min_sample_size = 100`, so `assemble_phyloseq()` returns NULL on both sides |
| different | 9 | SM_RefSeq_20250528 in every dataset: the intended fix (below) |

- Samples, taxa (except RefSeq, below) and `tax_table` are the same everywhere. Raw alpha diversity (Shannon, richness) is equal for every object other than RefSeq.
- The full GTDB rs220 gathers (`SM_gtdb_rs220`, P19 only) are no longer parsed, by decision.

**SM_RefSeq_20250528.** Same samples. Today's counts are 2.03–2.25× (median per dataset) the gather counts, which the merged objects now equal. The inflation isn't uniform: genomes are matched to 1, 2 or more lineage rows, and a few split across 2 species. So relative abundances shift too:

| Dataset | Taxa ref / new | Count ratio ref/new | Max Shannon diff | Max richness diff |
|---|---|---|---|---|
| AD_Skin | 338 / 338 | 2.05 | 0.33 | 1 |
| Bee | 25 / 25 | 2.25 | 0.15 | 0 |
| Maple | 731 / 729 | 2.03 | 0.13 | 2 |
| Moss | 396 / 396 | 2.22 | 0.74 | 1 |
| NAFLD | 338 / 338 | 2.05 | 0.37 | 0 |
| P19_Gut | 321 / 321 | 2.08 | 0.50 | 0 |
| P19_Saliva | 249 / 249 | 2.10 | 0.12 | 0 |
| PD | 874 / 874 | 2.06 | 0.52 | 1 |
| RA_Gut | 620 / 619 | 2.07 | 0.49 | 1 |

The Maple and RA_Gut taxa found only in the reference are species that reached `min_taxa_count` only through the doubled counts, or that got their counts from a genome split across 2 species.

**Fixed during the full run:**
- **AD_Skin `SM_gtdb-rs214-full/SRR20002431_gtdb-rs214-full_gather.csv` is truncated** at 48 KiB, mid-row 80 (an interrupted write, 2024-10-11). It's the only one of 6,812 gathers. sourmash can't load the partial row, so merge_tables.sh drops it from a temp copy. parse_SM() had put that row under a bogus sample "SRR2000", which assemble_phyloseq() dropped, so the objects match. The sample's rs214-full profile is still incomplete (79 matches kept): rerun that gather if it matters.
- **combine_mpa.py duplicates rows** when a taxon is a root in one input and comes after a prefix line in another (RA_Gut KB10: some per-sample files also have higher-rank lines). The 91 duplicates are exact copies. merge_tables.sh keeps the first, and stops if a taxon ever appears twice with different values.

**Noticed, not changed:**
- `data/PD/MOTUS/` holds 40 empty profiles named after the Bee samples (and 70 empty PD profiles). They're ignored on both sides.
- In 0_proc_data.R, `ps_raw.ls$Moss$MOTUS <- NULL` does nothing: the object is called `MOTUS3`.
- Maple has no sample data in 0_proc_data.R yet. It was compared with every sample of its tsv.
