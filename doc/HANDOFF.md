# Handoff between the HPC and Mac sessions

Read at the start of each session and update before stopping (see `CLAUDE.md`). Pending work items are in `doc/TODO.md`; this file says where things stand and what the last session did.

## Current state (keep this true; rewrite rather than append)

- **Branch:** `merged_tables`. Last pushed commit: `a8dfc1ad5`. 13 commits ahead of `main`; not merged yet.
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
- **Synced to the Mac:** done 2026-11-04. The Mac `Out/_Rdata/*.RDS` are the HPC files from 09:10, and all 9 datasets load.
- **No pre-fix baseline:** the old Mac RDS files were overwritten by the sync without a copy, and they were never in git. What's left: `Out/_Rdata/alpha_div_plot_data.tsv` and `alpha_div_wilcox.tsv` (2026-03-10), and possibly Time Machine (not mounted). The old vs new comparison is on hold until the user asks for it.
- **Figure scripts:**
  - `6.0`–`6.4` read only the synced files.
  - `6._FigX_TaxAssignment.R` and the DAA scripts (`5.x`) read files nobody regenerates anymore (`tax_assign_*`, `ps_taxranks`, `ps_filt`). DAA is out of scope.
  - The figure scripts are not in git (untracked on the HPC). The Mac copy is the reference.
- **Object names changed:** the sylph entries are `SYLPH_RefSeq_Seq/_Tax` and `SYLPH_gtdb_rs220_Seq/_Tax` (labels in `CCE_names`/`CCE_metadata` in `0_config.R`). The Moss MAGs object is `SM_gtdb-rs214-rep_MAGs`. Figure code that hard-codes old database names may need updating.

## Next steps

1. **[User, then either]** Find and process new public datasets. The candidate list is in `doc/ena_candidates/wgs_candidates_ranked.tsv` (regenerated 2026-11-04 by `scripts/src/ena_search/`). Then follow "Adding a dataset" in `doc/TODO.md`. The user will plan the next work after this.
2. **[Mac, on hold]** Compare the old and new results (alpha, Procrustes, PERMANOVA), focusing on sourmash RefSeq and MetaPhlAn, for the reviewers' response. The MetaPhlAn impact numbers measured so far are in `doc/TODO.md`.
3. **[Mac]** Rerun the figures (`6.0`–`6.4`). Decide how sylph appears in them (TODO: "Figures that include sylph need a rethink").
4. **[Either]** Merge `merged_tables` into `main` when the figures look right.

## Log (newest first; tag each entry [HPC] or [Mac])

### 2026-11-04 [HPC] (afternoon)
- **Mosquito (PRJNA1065965):** dropped, because of its paired/longitudinal design. Its partial download is gone; delete `data/Mosquito/` when convenient.
- **New dataset `Oilfield` (PRJNA1196024):**
  - 87 oil-contaminated coastal wetland soils: Shengli 44 (`S#-SL`), Dagang 43 (`D#-DG`).
  - NovaSeq 2×150, about 24.5 M pairs per sample, 282 GB.
  - Paper: Commun. Earth Environ. 2025, s43247-025-02062-1. The paper compares the two oilfields, so the grouping variable is `Oilfield`. Per-sample TPH and Fe values are on Figshare (paper's Data availability). Code: "on request".
  - **Download:** Slurm job 404290 on c-3t (cp0359), running `data/Oilfield/raw/download_ena.sh`: 8 samples in parallel from ENA FTP, with resume and md5 checks. Results go to `download_status.tsv` (174 files expected). The rate is about 2.5 MB/s in total, so it should take about 32 h. If the job stops early, resubmit it as is: finished files are skipped.
  - **Metadata:** `data/Oilfield/raw/metadata.tsv` (Sample = SRR, Alias, Site, Oilfield).
  - **Next:** kneaddata with the default GRCh38 index (soil, no host), then the "Adding a dataset" steps.
- **New dataset `Grapevine` (PRJNA1153634):**
  - Genome Biol 2026, 10.1186/s13059-025-03905-y. 78 runs from 39 vines: Health H/UH × Site A/D × Cultivar L/Z × Compartment Root/Soil. Each vine gives one root and one soil sample, so the compartments are paired.
  - All 78 runs are being downloaded (user's choice), 625 GB: job 404291, which starts after 404290 (`afterany`).
  - **Metadata:** `data/Grapevine/raw/metadata.tsv`, with Health, Site, Cultivar, Compartment, and Plant (the vine ID, for pairing).
  - **Root samples need a *Vitis vinifera* host index** before kneaddata. Check that R = root from the host fraction.
  - The download script is now generic: `data/<DS>/raw/download_ena.sh <filereport.tsv>`.
- **Downloads from NCBI's AWS copy, not ENA:** ENA FTP gave about 0.4 MB/s per stream. Three scripts in `data/<DS>/raw/`, chained as Slurm jobs:
  1. `download_sra.sh <filereport>`, on bio,c-3t (internet): `.sra` files into `sra/`, 8 at a time, checked against the server's size.
  2. `convert_sra.sh <filereport>`, an array with one run per task, on bio,c-3t,cpubase: `fasterq-dump` on local scratch, then pigz into `<run>_1/_2.fastq.gz` with ENA-style read names. Spots are checked against the ENA `read_count`. Writes `download_status.tsv`.
     - The first attempt ran 8 conversions in a single 32G job and was OOM-killed.
  3. `prep_kneaddata.sh <index>`: `samples_to_process.tsv` from the runs with status ok, then generates and submits the kneaddata array (60G, 24 threads).
  - Jobs (download → convert → prep):
    - Snake (3 test runs): 404307 → 404308 → 404309
    - Oilfield: 404310 → 404311 → 404312
    - Grapevine: 404313 (after 404310) → 404314 → 404315
  - Host indexes are built: `$FAST/host_genomes/VitisVinifera_index/VitisVinifera` (GCF_030704535.1) and `DrymarchonCorais_index/DrymarchonCorais` (GCA_043091225.1). Oilfield uses GRCh38.
  - Snake metadata: `data/Snake/raw/metadata.tsv` (Origin Wild/Captive from the BioSample location, all 80).
- **Grapevine:** all 78 runs converted (78 ok, 0 bad). Kneaddata array 405234 (78 tasks, Vitis index) submitted. Oilfield kneaddata is array 404625 (87 tasks, GRCh38).
- **Test plots (`scripts/tests/test_profiles.sh` + `test_composition_barplot.R`, outputs in `Out/Snake_test_3samples/[<DS>/]`):** Oilfield and Grapevine, 5 samples each, are being profiled (sylph 405312/405318, MetaPhlAn 405317/405319). Snake: kneaddata kept only 0.8–4.7% of read pairs (the rest is host). sylph detected 0% (wild zr21461_8), 28% and 50% of reads. `sylph-tax taxprof` crashed with 'Illegal instruction' on cp0404, so it was rerun on the login node. Watch for this node issue whenever sylph-tax runs inside Slurm. Barplot done: 92–96% host, 3% low quality; MetaPhlAn classifies 4–27% of the reads kept (Morganella, Providencia, Citrobacter, Yokenella…), sylph 0–50%. The wild sample is nearly empty. Whether the Snake dataset is usable is the user's call.
- **Candidates assessed, not downloaded:**
  - PRJNA1171820 (pond mesocosms): longitudinal, 2.25 TB. Rejected.
  - PRJNA1222911 (macaque Blastocystis): good; see the chat for details.
  - PRJNA1347217 (tomato Ralstonia): longitudinal, but one time point could be used. Its metadata is in the BioSample attributes.

### 2026-11-04 [Mac]
- Sync checked: the RDS files are the HPC 09:10 versions, with 9 datasets. No pre-fix copy was kept (see Current state).
- ENA search rerun. `find_wgs_projects.sh` and `rank_wgs_projects.R` now write to `doc/ena_candidates/` (untracked) instead of `Out/ena_candidates/`. The run found 710,018 runs (first public 2023+).
- Uncommitted Mac edits, not pushed yet:
  - `scripts/0_config.R`: Maple in `dataset_names`, shorter sylph labels and short names, DAA `tool_vars`/`DAA_metadata` commented out.
  - `scripts/2_alpha_compare.R`, `scripts/3_beta_compare.R`: Maple and the sylph objects added.

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
