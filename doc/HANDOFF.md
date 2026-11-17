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

**Done 2026-11-12** (job 408610, 26 min; log `logs/run_all_20261112-408610.log`). Checks passed:
- Moss has 140 samples, except MAGs 131 and the databases that already had fewer (MOTUS4 15, sylph 139/100). The Wilcoxon n stays at 131.
- PERMANOVA still has 221 rows.
- Oilfield, Grapevine and Horse have `SM_gtdb-rs220-rep`.

Classification rates for mOTUs, MetaPhlAn and sylph were rerun for all datasets. PD mOTUs 3 has 1 NA row left, `SRR19064900`, an empty 2025 file of a sample with no KneadData output, so not post-QC. `1b_beta_sensitivity.R` was resubmitted on the new `ps_rare` (job 408619, cpubase). It failed after the BC_VST PCoA: `0_config.R` now attaches `conflicted`, so `set_names` (purrr vs magrittr) is ambiguous. Lines 74–75 now call `purrr::set_names` (the only ambiguous name in the script; checked 2026-11-16). Resubmitted 2026-11-16 as job 408719, which completed in 34 min. The 6 objects `{pcoa_,procrustes_,permanova_}{BC_VST,rAitchison}` are written, and `sync_from_hpc.sh` already pulls them.

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

## Open decisions (HPC session, 2026-11-11/12)

- **Rarefaction depth: decided 2026-11-16, keep the minimum** (what job 408610 used, so no rerun). A sensitivity analysis with the outlier rule may come later. Measured 2026-11-16 on `ps_raw.ls` (rarefied tables only): the outlier rule would drop 219 sample × table cells, in 43 of 163 tables. Most are in Maple (88 cells, 32 samples; KB90 alone loses 29 of 88), AD_Skin (79 cells, 10 samples, 13 of 15 tables) and RA_Gut (26). Bee, Grapevine, NAFLD and both P19 datasets lose none. Rarefaction is per dataset × database on classified counts, at the minimum sample. Options:
  - keep the minimum;
  - use the 5th percentile, which drops about 5 % per table: 1–4 samples in most datasets, 8 in RA_Gut, 35 in PD, 5–37 in Horse;
  - use an outlier rule: drop a sample from one table only below 10 % of that table's median.
  - Decided: no post-QC depth floor; sparse pairs are kept and footnoted. The minimum means rerunning only the 4 datasets above; the other options mean editing `0_proc_data.R` and rerunning all 12.
- **Sensitivity objects** `pcoa_/procrustes_/permanova_{BC_VST,rAitchison}` (from `scripts/1b_beta_sensitivity.R`, untracked) date from the 2026-11-09 `ps_rare`. Rerun the script after `0_run_all.R`.
- **Supplementary sample table** (planned): per sample, post-QC pairs and reads, metadata/group, the reason it was dropped, and presence per database after rarefaction.

## Elsewhere

- **Resource benchmark, run on the HPC** (not Fir), on cp0359 only, so on uniform hardware. 55 samples (no Horse); job ids in `Out/benchmark/jobs.tsv`.
  - Done and checked: bases, sylph and sourmash sketches, SYref, SY220, MPA22/23/25, MOTU3.
  - **KB and KG (408048/9): 35 of 55 samples ran on 0 reads** (rows 1–10, 16–40), though they show as COMPLETED. Their `tools.tsv` paths started with `/net/nfs-ip34/`, the same filesystem as `/net/nfs-bio/`, but `kraken_local.sh` only makes `$ILAFORES` visible in the container. `tools.tsv` was fixed on 2026-11-12 (backup: `tools.tsv.bak_nfs-ip34`). Rerun on c-3t with KB 408612 and KG 408613. Those also gave empty kreports (moved to `Out/benchmark/_empty_kreports_20261112/`), so the rows were rerun again with KB 408620 and KG 408621. Checked 2026-11-16: all 55 kreports in each of the 6 KB*/ folders have reads. For these reruns, `benchmark_kraken.sh` takes an optional 7th argument, a row list.
  - MOTU4 tasks 30 and 40 ran out of memory at 100G (bwa killed, then a ZeroDivisionError). Rerun at 200G with 408611: row 40 completed (100G peak), but row 30 ran out of memory again. That row is PD SRR19064482, the largest sample (2 × ~56 GB fastq). Resubmitted 2026-11-16 at 400G as 408720 (array=30): completed, peak 277G.
  - `10_benchmark_plot.R` now reads every job of a label and keeps, per label × row, the last COMPLETED job. Without this, a rerun would hide the rest of the label, and the 0-read steps would enter the figure.
  - All benchmark jobs stay on c-3t (cp0359) so the hardware is the same for every tool (decided 2026-11-12). Other heavy work goes to cpubase.
  - **Benchmark complete (2026-11-16):** all 17 tools have 55/55 samples. `sacct.tsv` was regenerated and `10_benchmark_plot.R` run on the HPC, which wrote `resource_use.tsv`, `summary.tsv`, `benchmark_ridges.pdf` and `benchmark_scaling.pdf` in `Out/benchmark/`.
  - `doc/HANDOFF_FIR.md` describes the original Fir plan.
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
