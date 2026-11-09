You're picking up work on the MethodsComparison repo. It compares taxonomic profilers (mOTUs 3/4, MetaPhlAn 4, Kraken2/Bracken, sourmash, sylph) on several metagenomic datasets. The earlier sessions ran on my Mac, so I'm summarising where things stand. Read the files I point to before acting. Treat any Mac path (`/Users/jorondo/...`, `~/ip34`) mentioned anywhere as stale.

## Setup

- **This clone:** the HPC clone at `$MC` = `/net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison` (GitHub `jorondo1/MethodsComparison`, branch `main`). Profilers run here through SLURM. All the commands are documented in `scripts/src/MethodsComparison_Pipeline.sh`. It shows how each tool and environment is invoked: motus venv, `ml apptainer` for Kraken, sourmash, the sylph conda env under `$ILAFORES/programs/conda_envs/`.
- **The Mac clone:** this is where R runs (`scripts/0_proc_data.R` and onward). It gets the result files via rsync, not git: `data/` and `Out/` are now gitignored.
- **Git history:** it still contains old data. The clones were synced with `git fetch` + `git reset origin/main` (a mixed reset). Never run `git reset --hard` or `git clean`, which would destroy untracked results in `data/`. Don't push without asking me.
- **Compute:** run only light commands on the login node. Use `sbatch` or an interactive allocation for anything heavy or long. The old compute nodes can't see `~/miniconda3`, so put environments under `$ILAFORES`.
- **Never rerun the profilers** for this task. Work only from the existing per-sample outputs.

## Already done (don't redo or "fix")

- **sylph 1.0.0:**
  - Scripts: `scripts/src/sylph_SLURM.sh` profiles each sample into `data/<ds>/SYLPH_<db>/per_sample/`. `scripts/utils/sylph_merge.sh [-f] <ds> <db> <tsv>` writes `counts.tsv`, `taxonomy.tsv`, `merged_species.tsv` and `classification_rate.tsv`.
  - Counts are divided by the classification rate on purpose. `sylph profile --estimate-read-counts` scales counts by the rate twice, so they come out as rate² × reads. This was checked empirically.
- **MetaPhlAn normalisation:**
  - `parse_MPA(MPA_data = TRUE)` in `jorondo1/misc_scripts/community_functions.R` now normalises per sample (commit b86aac6). `0_proc_data.R` sources it from GitHub main.
  - The old version divided by the sum over all samples. The results haven't been regenerated yet (I know; that's on my list).
- **Classification rate fixes** in the pipeline: mOTUs uses `awk -F'\t'` and rate = assigned/total. MetaPhlAn uses unclassified = 100 − Σ species rows; the current profiles have no UNCLASSIFIED line and no higher ranks.

## The task: merged-tables side project (on a new branch)

**Goal:** stop reading thousands of per-sample files in R. For each tool, merge the per-sample profiles on the HPC with that tool's own merge utility, then adapt the R parsing to read the merged tables. The resulting phyloseq objects must be the same as the current ones.

### Versions that produced the current per-sample files

Use these exact versions for the merge utilities:

| Tool | Version | Merge utility | Note |
|---|---|---|---|
| mOTUs | 3.1.0 and 4.0.4 | `motus merge` | |
| MetaPhlAn | 4.2.4 | `merge_metaphlan_tables.py` | |
| sourmash | 4.8.11 | `sourmash tax metagenome` | |
| Kraken2 | 2.1.2 | KrakenTools 1.2 `combine_mpa.py` | Bracken 2.8 also has `combine_bracken_outputs.py`, but it carries no lineage, so prefer `combine_mpa.py`. |

### How R reads each tool today (this defines "the same")

Everything is in `scripts/0_proc_data.R` (`meta_parsing()`) plus `misc_scripts/community_functions.R`: https://github.com/jorondo1/misc_scripts/blob/main/community_functions.R. Read both.

**mOTUs**
- Files: `data/<ds>/MOTUS/*_profile.txt` (v3; columns mOTU, Taxonomy, NCBI, Abundance) and `data/<ds>/MOTUS4/*_profile.txt` (v4; columns mOTU, Taxonomy, Abundance).
- Lines starting with `#` are dropped, and so is mOTUs 4's bare `Taxonomy` header line. Rows containing `unassigned` are dropped.
- Only species rows are kept (`s__` without `t__`). The species key is `Taxonomy` pasted with the mOTU ID, so that several "incertae sedis" mOTUs stay distinct. Values are counts.

**sourmash**
- Files: `data/<ds>/SM_<db>/*_gather.csv`.
- Count per genome = `unique_intersect_bp / scaled * average_abund`, rounded.
- The genome ID is `match_name` (or `name` in older outputs): first token, `.fa.sig` stripped, then the part between the first `_` and the first `.`. Plasmid matches are dropped.
- Taxonomy comes from `SM_<db>/SM_<db>_lineages.csv`. That file is a headerless subset of the DB lineage CSV under `$ILAFORES/ref_dbs/sourmash_db/`, so `sourmash tax` should use the full DB lineage files instead. It's parsed by `parse_GTDB_lineages` or `parse_genbank_lineages`, then `species_glom` sums the genomes into species.

**Kraken/Bracken**
- Files: `data/<ds>/KB*/*/*_bracken/*_bracken_S.MPA.TXT` (`kreport2mpa` output of the species-level Bracken reports).
- Columns are Taxonomy and Abundance. The sample name comes from the file name.

**MetaPhlAn**
- Files: `data/<ds>/MPA_*/*/*_profile.txt` (columns Taxonomy, NCBI, Abundance, Void).
- Only species rows are kept (`s__` without `t__`), then each sample is renormalised to 100.

**sylph**
- Already merged: `counts.tsv` + `taxonomy.tsv`, read by `parse_sylph()`.

**All tools**
- Everything goes through `mgx.tools::assemble_phyloseq(samData, seqdepth_plot = FALSE, onlySpecies = TRUE)` (https://github.com/jorondo1/mgx.tools). MetaPhlAn also gets `min_sample_size = 0, min_taxa_count = 0`.
- Special cases: Moss drops some databases, and Moss's `SM_gtdb-rs214-rep_MAGs` comes from an external RDS. These stay as they are.

### Open questions to settle first

Test on 2–3 samples of one dataset before running at scale.

1. **MetaPhlAn:** what `merge_metaphlan_tables.py` (4.2.4) names the sample columns: does it strip `_profile`? What header and comment lines does it write? Does it cope with profiles that lack the UNCLASSIFIED line and higher ranks?
2. **mOTUs:** for both 3.1.0 and 4.0.4, does `motus merge` keep the mOTU IDs and the full `k__…|s__…` taxonomy the profiles were made with? Check the profiling flags in `scripts/src/motus3_SLURM.sh` and `motus4_SLURM.sh`. Does the merge need the database installed? Are the values counts?
3. **Kraken:** where does `combine_mpa.py` get the sample names (a header line in the input, or file order)? Are missing taxa filled with 0?
4. **sourmash:** `tax metagenome -F csv_summary` with the gather CSVs and the full DB lineage CSV.
   - Does `f_weighted_at_rank × total_weighted_hashes` at the species rank reproduce the summed `unique_intersect_bp / scaled * average_abund`?
   - How are matches without a lineage (plasmids) reported?
   - Does identifier matching work for each DB (gtdb-rs214 reps/full, gtdb-rs220 reps, genbank-2022.03, RefSeq_20250528)?
   - How are NCBI lineages written (with or without rank prefixes)?
   - `tax metagenome` also reports an unclassified fraction, which might replace the sourmash classification-rate block.

**Recommended order:** sourmash first (it has the most R gymnastics to remove), then the MetaPhlAn-style tools (MetaPhlAn, Kraken, mOTUs).

### Validation

R and phyloseq may not be set up on the HPC; check with `module avail` and whether `mgx.tools` is installed. If they aren't, produce the merged tables and the R changes here, and I'll run the comparison on my Mac. The comparison, for every dataset × database:

- **Must match:** same samples, same taxa and the same `tax_table`. `otu_table` values must be equal; report and explain any difference, e.g. rounding in the sourmash counts.
- **Diversity:** alpha diversity on the raw objects should be identical.
- **Rarefied objects** may differ: rarefaction depends on taxa order, even with seed 1234. I accept that.
- Write the comparison as a script, so it can be rerun.

### Deliverables, on branch `merged_tables` (create it from `main`)

- Merge commands in `scripts/src/MethodsComparison_Pipeline.sh`, as one block per tool looping over datasets and DBs (like the sylph merge loop), or as small scripts like `sylph_merge.sh`. Write the merged tables next to the per-sample files, e.g. `data/<ds>/<tool dir>/merged_*.tsv`, and propose a naming scheme. Never delete per-sample files.
- New R readers, roughly two:
  - one for MPA-style wide tables with a lineage string;
  - one for sourmash `csv_summary`.
  
  `meta_parsing()` would then become a loop over tool → files → reader. Keep the `assemble_phyloseq` calls.
- The comparison script and a short report of the results.
- My Mac rsync script (`scripts/sync_from_hpc.sh`, untracked, Mac only) will need its include list switched to the merged files. Tell me which patterns.

Commit on the branch as you go, but ask me before pushing. We'll merge the branch into main once everything checks out.

## Not part of this task (don't start without asking)

- Regenerating the results after the MetaPhlAn fix.
- Rerunning sylph and the classification rates.
- The sylph figures.
- The git history purge: removing `data/` and `Out/` from history before the next Zenodo snapshot. It's deferred, and the steps live on my Mac.

Start by reading `scripts/0_proc_data.R`, `community_functions.R`, the parts of the pipeline that run each tool, and one example per-sample output per tool (pick a small dataset). Then tell me your plan for the format checks before you start installing or running anything.
