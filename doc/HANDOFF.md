# Handoff between the HPC and Mac sessions

Read at the start of each session and update before stopping (see `CLAUDE.md`). Pending work items are in `doc/TODO.md`; this file says where things stand and what the last session did.

## Current state (keep this true; rewrite rather than append)

- **Branch:** `merged_tables`. Last pushed commit: `01e908d2e`. 12 commits ahead of `main`; not merged yet.
- **Pipeline (all on the HPC):**
  - Merged tables per tool, from `scripts/src/merge_tables.sh` and `sylph_merge.sh`.
  - `scripts/0_run_all.R` runs, per dataset, `0_proc_data.R` (phyloseq objects + rarefaction) and then `1_diversity_metrics.R`. Results are cached in `Out/_Rdata/cache/<DS>/`, then assembled into `Out/_Rdata/*.RDS`.
  - The "Adding a dataset" steps are in `doc/TODO.md`.
- **Results in `Out/_Rdata/` on the HPC:** regenerated on 2026-11-04 (09:10) for 9 datasets: P19_Saliva, P19_Gut, Moss, NAFLD, AD_Skin, RA_Gut, Bee, PD, Maple. They include:
  - the MetaPhlAn per-sample normalisation fix;
  - the sourmash RefSeq double-counting fix (`doc/sourmash_refseq_double_counting.txt`);
  - Bacteria + Archaea only (GTDB tables kept whole);
  - the new sylph objects `SYLPH_<db>_Seq` (rarefied counts) and `SYLPH_<db>_Tax` (relative abundance rescaled to 100, not rarefied), for RefSeq and gtdb_rs220;
  - Maple (grouping variable `Provenance`: 88 samples, 66 tested) and Moss sylph;
  - seeded Procrustes and PERMANOVA.
- **Synced to the Mac:** NOT YET (as of 2026-11-04).
- **Mac `Out/_Rdata/`:** still the old results from before the fixes. Copy it aside first (e.g. to `Out/_Rdata_preFix_<date>/`); it's the baseline for measuring how the fixes change the results.
- **Figure scripts:**
  - `6.0`–`6.4` read only the synced files.
  - `6._FigX_TaxAssignment.R` and the DAA scripts (`5.x`) read files nobody regenerates anymore (`tax_assign_*`, `ps_taxranks`, `ps_filt`). DAA is out of scope.
  - The figure scripts are not in git (untracked on the HPC). The Mac copy is the reference.
- **Object names changed:** the sylph entries are `SYLPH_RefSeq_Seq/_Tax` and `SYLPH_gtdb_rs220_Seq/_Tax` (labels in `CCE_names`/`CCE_metadata` in `0_config.R`). The Moss MAGs object is `SM_gtdb-rs214-rep_MAGs`. Figure code that hard-codes old database names may need updating.

## Next steps

1. **[Mac]**
   1. `git pull` (branch `merged_tables`).
   2. Copy `Out/_Rdata/` aside.
   3. `bash scripts/utils/sync_from_hpc.sh <user@host>`.
   4. Check that the RDS files load and have 9 datasets.
2. **[Mac]** Compare the old and new results (alpha, Procrustes, PERMANOVA), focusing on sourmash RefSeq and MetaPhlAn, for the reviewers' response. The MetaPhlAn impact numbers measured so far are in `doc/TODO.md`.
3. **[Mac]** Rerun the figures (`6.0`–`6.4`). Decide how sylph appears in them (TODO: "Figures that include sylph need a rethink").
4. **[Either]** Merge `merged_tables` into `main` when the figures look right.

## Log (newest first; tag each entry [HPC] or [Mac])

### 2026-11-04 [HPC]
- **sylph:**
  - `sylph_merge.sh` now also writes `merged_species_tax.tsv` (taxonomic abundance), plus a `--tax-only` mode.
  - Ran it for all datasets.
  - Profiled Moss with sylph, RefSeq + GTDB (first sylph run for Moss).
  - Moss `datasets.tsv` now points to `data/Moss/raw/preprocessed_reads.sample.tsv`, with the corrected read paths.
- **`scripts/utils/parse_merged.R`:**
  - `prokaryotes_only()`: Bacteria + Archaea; Tax/MPA rescaled to 100.
  - `parse_sylph(abundance =)`.
  - `dataset_fixes()`: Moss MAGs object; Maple renames.
  - sylph fungi_refseq and glob_r232 are excluded (side tests only).
- **Per-dataset restructure:**
  - `0_proc_data.R` and `1_diversity_metrics.R` now take a dataset argument.
  - New master script `scripts/0_run_all.R`: checks the setup of new datasets, caches, assembles.
  - Verified: the assembled outputs are `identical()` to the previous monolithic run. The only differences are the p-values, from the new seeds.
- **Diversity:** samples without a group are dropped only from the Wilcoxon tests and PERMANOVA. `set.seed(1234)` added before `adonis2` and before Procrustes.
- **Archived** (not in the repo): the old-vs-new merged-table comparison scripts (`archive/merged_tables_comparison/`).
- **Docs:**
  - `doc/sourmash_refseq_double_counting.txt` (for the reviewers);
  - `doc/TODO.md` methods notes (the test sample sets differ from the diversity sample sets);
  - the "Adding a dataset" checklist.
- **`scripts/utils/sync_from_hpc.sh`:** now also syncs the sylph tax tables and the assembled RDS files.
- **Stale leftovers** on the HPC in `Out/_Rdata/`: `ps_filt`, `ps_radEmu`, `ps_taxranks`, `ps_rare_taxranks` (June/July 2025). They are not synced; deleting them is the user's call.
