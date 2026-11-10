# Handoff

Read at the start of a session and keep it current. Pending items are in `doc/TODO.md`. Older session logs are in git history (`git log -p -- doc/HANDOFF.md`).

## Where things stand (2026-11-10)

- **HPC work is finished.** Profiling, merged tables, classification rates and `scripts/0_run_all.R` are complete for all 12 datasets:
  - P19_Saliva, P19_Gut, Moss, NAFLD, AD_Skin, RA_Gut, Bee, PD, Maple;
  - the new Oilfield, Grapevine and Horse.
  - From here the work is on the Mac, from the synced outputs. The pipeline runs only on the HPC; rerunning it is not planned.
- **Results:** `Out/_Rdata/*.RDS`, assembled on 2026-11-09 at 22:26. Per-dataset caches are in `Out/_Rdata/cache/<DS>/`.
  - **PD** was redone with the re-merged mOTUs 3 table (724 samples). SRR19064900 has no preprocessed reads in any tool.
  - **Grouping variables:**
    - Oilfield: `Oilfield`.
    - Grapevine: `Health`, also has `Compartment` and `Plant` for pairing.
    - Horse: `island_side`, also has `age_class_rev`.
  - **Samples after rarefaction:**
    - Oilfield: 87, except mOTUs 3 (81) and sylph RefSeq (85).
    - Grapevine: 78 everywhere.
    - Horse: 744 for most databases. RefSeq covers horse gut poorly: mOTUs 3 587, mOTUs 4 682, KB90 525, sourmash RefSeq 182, sylph RefSeq 140/144.
- **Classification rates:** `Out/classification_rates/*.tsv`, all 12 datasets. Horse has no `reads` rows, because its KneadData logs stayed on Fir.
  - New `kraken_species_classification_rate.tsv`: Kraken's own species-rank share, read from the Bracken summaries. Use it instead of `kraken`/`bracken` for species-level identifiability.
  - mOTUs is dropped from the rate figure: its denominator is marker-gene reads only. TODO "Methods notes" has the wording item.
- **Object names:** the sylph entries are `SYLPH_<db>_Seq` (rarefied counts) and `SYLPH_<db>_Tax` (relative abundance, not rarefied). The Moss MAGs object is `SM_gtdb-rs214-rep_MAGs`.
- **Git:** branch `benchmark` (from `main`).
  - Staged: `scripts/benchmark/`, `scripts/10_benchmark_plot.R`, `doc/HANDOFF_FIR.md`.
  - Unstaged: the R setup for the 3 new datasets (`0_config.R`, `0_run_all.R`, `utils/sample_data.R`), `src/classification_rates.sh`, `utils/sync_from_hpc.sh` and docs.
  - Untracked: `src/kraken_node_SLURM.sh`.
  - Nothing is committed. The figure scripts (`6.x`, `8_*`, `9.x`) are untracked on the HPC; the Mac copy is the reference.

## Sync to the Mac

From the Mac repo root:
```bash
bash scripts/utils/sync_from_hpc.sh <user@host>      # add --dry-run to preview
```
It pulls only what the figure scripts read:
- `Out/classification_rates/*.tsv`;
- `Out/_Rdata/{ps_raw.ls,ps_rare.ls,alpha_div,pcoa_noVST.ls,pairwise_dist,procrustes,permanova.ds,dataset_names_n_lookup_alphadiv}.RDS`.

`-u` keeps local files that are newer. `scripts/src/datasets.tsv` comes with git. Nothing else on the HPC is needed on the Mac.

## Next steps (Mac)

1. Sync, then rerun the figures `6.0`–`6.4` with the 12 datasets. Add display names or ordering for Oilfield, Grapevine and Horse where scripts hard-code dataset lists.
2. `8_Identifiability.R`: add `'Horse'` to `dataset_order` (line 26) and to `ds_names`. Otherwise `stopifnot(setequal(...))` fails.
3. `8_Classification_rates.R`: remove mOTUs, and switch Kraken to `kraken_species_classification_rate.tsv`.
4. Decide how sylph appears in the figures (TODO).
5. Stale, not regenerated: `6._FigX_TaxAssignment.R`, `9.2` (`Rarefaction.rds`), DAA `5.x`. Out of scope unless the user says otherwise.

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
