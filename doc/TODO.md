# MethodsComparison TODO (untracked: doc/ is in .gitignore)

Started 2026-11-03.

## Results to regenerate

- [ ] **MetaPhlAn per-sample normalisation fixed** (misc_scripts `b86aac6`, pushed 2026-11-03).
  - **The bug:** `parse_MPA(MPA_data = TRUE)` divided by the sum over all samples. Each sample therefore summed to its own classified fraction instead of 100.
  - **To do:** rerun the figures (objects and diversity regenerated 2026-11-04 with `scripts/0_run_all.R`).
  - **Impact measured** on the saved RDS (MetaPhlAn objects only):
    - Alpha diversity is unchanged: every index normalises per sample.
    - Bray-Curtis old vs new correlation: ≥ 0.997 for P19 Saliva/Gut, NAFLD, RA_Gut and PD; 0.979 for Moss; ~0.93 for AD_Skin.
    - PERMANOVA: R² changes in the 3rd digit, and no p-value crosses 0.05.
    - Procrustes: 540 MetaPhlAn pairs; 52 change by > 0.02 and 8 by > 0.05.
      - Moss, MPA_db2025 vs mOTUs4: 0.65 → 0.52.
      - AD_Skin, MPA vs KB90_GTDB: 0.85 → 0.91. Most other AD_Skin MPA pairs rise by 0.02–0.04; MPA vs KB10_GTDB falls from 0.95 to 0.91.
  - [ ] Check any figure or manuscript text about AD_Skin or Moss procrustes agreement involving MetaPhlAn.
- [ ] Classification rates: rerun the mOTUs block (fixed: `awk -F'\t'` and rate = assigned/total), the MetaPhlAn block (100 − Σ species) and the sylph block in `MethodsComparison_Pipeline.sh`.

- [ ] **sourmash RefSeq counts are doubled** (found 2026-11-03, see `doc/sourmash_refseq_double_counting.txt`). R joins gathers to `SM_RefSeq_20250528_lineages.csv` on the bare accession number, which matches the GCA/GCF twins and versions, so `left_join` duplicates most genomes (+111–116 % on P19_Saliva GQ10–12). Counts, rarefaction depth and count-based alpha are affected. Fixed on the `merged_tables` branch (RefSeq lineages from the sylph-tax table). The relative abundances shift too (max Shannon difference 0.1–0.7 per dataset), so the sourmash RefSeq results need regenerating.

## Methods notes (for the manuscript)

- [ ] **The sample sets differ between the diversity comparisons and the statistical tests.**
  - The diversity comparisons use every sample of a dataset: alpha diversity, PCoA, pairwise distances and Procrustes.
  - The Wilcoxon tests and PERMANOVA use only the samples that have a value for the grouping variable (`grouping_variable` in `0_config.R`).
  - Maple: 88 samples in the diversity comparisons, but only the 66 with a Provenance (33 M, 33 S) in the tests.
  - The `n` in the figure labels (`dataset_names_n_lookup_alphadiv.RDS`) comes from the Wilcoxon tests, so for Maple it's the tested count.
  - State this in the methods.
- Out of scope for this project: rarefaction curves (`9.1`/`9.2`) and the DAA scripts (`5.0`–`5.4`).

## Adding a dataset

`scripts/0_run_all.R` builds the phyloseq objects and the diversity metrics one dataset at a time. Results are cached per dataset in `Out/_Rdata/cache/<DS>/`, and the script then assembles the usual `Out/_Rdata/*.RDS` files. Only datasets without a complete cache are computed; list datasets in `rerun` to recompute them anyway.

1. HPC:
   - run the profilers;
   - add a line to `scripts/src/datasets.tsv`;
   - run `bash scripts/src/merge_tables.sh <DS>`, plus `scripts/src/sylph_merge.sh <DS> <db> <sample tsv>` for each sylph database.
2. `scripts/utils/sample_data.R`: add the sample data to `load_sample_data()`'s list.
3. `scripts/0_config.R`: add `grouping_variable` and `dataset_names` entries.
4. Add `<DS>` to `datasets` in `scripts/0_run_all.R` and run it on the HPC (`Rscript scripts/0_run_all.R`).
   - Before parsing, it checks steps 1–3 and reports anything missing, with the fix.
5. Sync (`scripts/utils/sync_from_hpc.sh`), then make the figures.

## sylph

- [ ] Rerun sylph 1.0.0 (env `$ILAFORES/programs/conda_envs/sylph`) for all datasets and DBs. First delete the old `per_sample/` dirs and the sketches made with 0.9.0.
- [ ] Merge with `sylph_merge.sh [-f]`.
- [ ] One-off: build `merged_species.tsv` from the existing `merged_all_ranks.tsv`, then delete `merged_all_ranks.tsv`.
- [ ] Check with the sylph authors that dividing the counts by the classification rate is correct (`doc/sylph_methods_note.txt`).
- [ ] Figures that include sylph: they need a rethink (only the parsing and config are done).

## Code issues noticed

- [x] `1_diversity_metrics.R` sourcing `myFunctions.R` / `0_Config.R`: fixed. It runs on the HPC.
- [x] `sync_from_hpc.sh`: now tracked, as `scripts/utils/sync_from_hpc.sh`.

## Merged tables (done, branch `merged_tables`)

All tools are read from the merged tables (`scripts/src/merge_tables.sh`, `scripts/utils/parse_merged.R`). The objects were checked against the per-sample parsing before the comparison scripts were archived. They were identical, apart from MetaPhlAn (float precision) and sourmash RefSeq (the double-counting fix).

## Before the next Zenodo snapshot

- [ ] Purge data/ and Out/ from the git history: `doc/git_history_purge_steps.txt`.
