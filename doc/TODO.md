# MethodsComparison TODO (untracked: doc/ is in .gitignore)

Started 2026-11-03.

## Results to regenerate

- [ ] **MetaPhlAn per-sample normalisation fixed** (misc_scripts `b86aac6`, pushed 2026-11-03).
  - **The bug:** `parse_MPA(MPA_data = TRUE)` divided by the sum over all samples. Each sample therefore summed to its own classified fraction instead of 100.
  - **To do:** rerun `0_proc_data.R`, then `1_diversity_metrics.R`, then the figures. Other tools' objects shouldn't change.
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

## sylph

- [ ] Rerun sylph 1.0.0 (env `$ILAFORES/programs/conda_envs/sylph`) for all datasets and DBs. First delete the old `per_sample/` dirs and the sketches made with 0.9.0.
- [ ] Merge with `sylph_merge.sh [-f]`.
- [ ] One-off: build `merged_species.tsv` from the existing `merged_all_ranks.tsv`, then delete `merged_all_ranks.tsv`.
- [ ] Check with the sylph authors that dividing the counts by the classification rate is correct (`doc/sylph_methods_note.txt`).
- [ ] Figures that include sylph: they need a rethink (only the parsing and config are done).

## Code issues noticed

- [ ] `1_diversity_metrics.R` sources `scripts/myFunctions.R`, which no longer exists (the functions moved to mgx.tools?). It also sources `scripts/0_Config.R` while the file is `0_config.R`. That only works on the Mac's case-insensitive disk and will fail on Linux.
- [ ] `sync_from_hpc.sh` (untracked): commit it, or keep it local.

## Side project: merged tables (branch, not started)

Merge each tool's per-sample profiles on the HPC with its own utility, read the merged tables in R, and check that the phyloseq objects come out identical (rarefaction differences are accepted).

- **Versions that produced the current tables:** mOTUs 3.1.0 and 4.0.4, MetaPhlAn 4.2.4, sourmash 4.8.11, Kraken 2.1.2, KrakenTools 1.2, Bracken 2.8.
- **Utilities:**
  - `merge_metaphlan_tables.py`
  - `motus merge`
  - KrakenTools `combine_mpa.py`
  - `sourmash tax metagenome -F csv_summary`
  - sylph: already done.
- **Start with sourmash:** it's where most of the R parsing work is.
- **Format questions to settle first:**
  - sample naming in the MetaPhlAn and Kraken merged tables;
  - whether `motus merge` keeps the mOTU IDs, for both mOTUs versions;
  - the format of NCBI lineages in sourmash `tax` output.

## Before the next Zenodo snapshot

- [ ] Purge data/ and Out/ from the git history: `doc/git_history_purge_steps.txt`.
