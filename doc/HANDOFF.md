# Handoff

Read at the start of a session and keep it current. Pending items are in `doc/TODO.md`. Older session logs are in git history (`git log -p -- doc/HANDOFF.md`).

## Where things stand (2026-11-10, end of a Mac session)

- **Results:** `Out/_Rdata/*.RDS` from 2026-11-09 (all 12 datasets), synced to the Mac. Classification rates for all 12, including Horse's KneadData `reads` rows.
- **Mac work, not committed yet** (the user commits and pushes before the HPC session):
  - `0_config.R`: dataset list and order (`these_datasets <- names(dataset_names)`) and shared helpers for the `6.x` scripts (`read_rates()`, `samples_qc()`, `samples_grouped()`, `dataset_names_n_qc()`, `rate_databases`). `0_run_all.R` sources it too, which is harmless: function definitions only.
  - Figure scripts `2_`–`5_` cover the 12 datasets. New `6.1_Identifiability.R` (ridges, replaces `8_*`) and `6.2_dataset_summary.R` (profiled-samples heatmap, `Out/dataset_summary/dataset_table.html`).
  - `utils/parse_merged.R`: `SM_gtdb-rs220-reps` → `SM_gtdb-rs220-rep` for every dataset (was Maple only), and the new sample rule (below). New tracked file `scripts/src/dropped_samples.tsv`.
  - `src/count_reads.sh` and `src/compile_sample_reads.jl` untracked and archived; `Out/stats/` is no longer used.
  - To stage: `git add -u && git add scripts/6.1_Identifiability.R scripts/6.2_dataset_summary.R scripts/src/dropped_samples.tsv && git add -f doc/HANDOFF.md doc/TODO.md`.

## Task for the HPC session: rerun `0_run_all.R` for Moss, Oilfield, Grapevine, Horse

The code is ready (Mac, 2026-11-10); the HPC only runs it. `scripts/0_run_all.R` already has:
```r
rerun <- c('Moss', 'Oilfield', 'Grapevine', 'Horse')
```
Then run it through sbatch, and sync to the Mac (`scripts/utils/sync_from_hpc.sh`).

**Why these 4:**
- **Oilfield, Grapevine, Horse:** sourmash GTDB 220 was named `SM_gtdb-rs220-reps`, so the figures left it out. `parse_merged.R` now renames it for every dataset.
- **Moss:** gets its 9 samples without metadata back (140 instead of 131). They're in the diversity comparisons and out of the tests (NA `Compartment`).

**New sample rule (decided with the user):**
- Every post-QC sample enters the diversity comparisons. Samples without metadata get NA groups and are only left out of the Wilcoxon tests and PERMANOVA.
- The exception is `scripts/src/dropped_samples.tsv`, which removes samples from everything:
  - P19_Gut: controls GQ27–29 and GQ54–56, and the mock community MSA_3001;
  - PD: the 5 SP/SC runs;
  - AD_Skin: the 78 runs without study metadata.
- It's applied in `meta_parsing()` (`utils/parse_merged.R`) and in `samples_qc()` (`0_config.R`, the `6.x` counts).
- Tested on the Mac for Moss, P19_Gut, AD_Skin and PD: only Moss changes. The others keep their 21 / 70 / 719 samples, so they need no rerun.

**Checks after the run:**
- Moss: 140 samples in most objects. `SM_gtdb-rs214-rep_MAGs` stays at 131: it's the external `mossMAGs.RDS` object.
- The Wilcoxon/PERMANOVA n are unchanged: Moss 131 in `dataset_names_n_lookup_alphadiv.RDS`.
- Oilfield, Grapevine and Horse have `SM_gtdb-rs220-rep`.

**Then, on the Mac:** rerun the figures `2_`–`6.2`. `6.2` writes `Out/dataset_summary/dataset_table.html` (2,216 samples post-QC, 2,185 with metadata).

## Sync to the Mac

From the Mac repo root:
```bash
bash scripts/utils/sync_from_hpc.sh <user@host>      # add --dry-run to preview
```
It pulls only what the figure scripts read:
- `Out/classification_rates/*.tsv`;
- `Out/_Rdata/{ps_raw.ls,ps_rare.ls,alpha_div,pcoa_noVST.ls,pairwise_dist,procrustes,permanova.ds,dataset_names_n_lookup_alphadiv}.RDS`.

`-u` keeps local files that are newer. `scripts/src/datasets.tsv` comes with git. Nothing else on the HPC is needed on the Mac.

## Next steps (Mac), after the rerun

- Render `2_`–`5_`, `6.1`, `6.2` with the new objects. `3_beta_compare.R` line 38 filters on `these_datasets` (from the config).
- Still open: how sylph appears in the figures (TODO). Stale, out of scope: `6._FigX_TaxAssignment.R`, `9.2`, DAA `5.x`.

## Elsewhere

- **Resource benchmark on Fir:** see `doc/HANDOFF_FIR.md`. The reads are staged as symlinks in `Out/benchmark/reads/` (55 samples; 5 Horse samples are added on Fir). The results come back as tsv files to `Out/benchmark/` and are plotted with `scripts/10_benchmark_plot.R`.
- **HPC cleanup (the user's call):**
  - stale `Out/_Rdata/{ps_filt,ps_radEmu,ps_taxranks,ps_rare_taxranks}` (2025);
  - `data/Grapevine/_partial_*`;
  - `data/Mosquito/`.
  - Jobs 405784_27–30 are stuck COMPLETING on cp0419; only admins can clear them.

## Rules

- No commits or pushes unless asked. Never `git reset --hard` or `git clean`. Sync clones with fetch + mixed reset.
- Don't rerun profilers unless asked. Never delete per-sample profiler outputs.
- Never replace a script that running jobs read (NFS: they stop silently with exit 0).
- HPC interactive node: at most 2 cores unless the user grants more. Heavy work goes through sbatch.
- Out of scope: GTDB 214, GenBank, full GTDB rs220, Olive.
