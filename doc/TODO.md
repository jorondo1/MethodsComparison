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
- [ ] Classification rates: rerun the mOTUs block (fixed: `awk -F'\t'` and rate = assigned/total), the MetaPhlAn block (100 − Σ species) and the sylph block: `rate_motus`, `rate_mpa`, `rate_sylph` (`scripts/src/classification_rates.sh`, called from `pipeline_2_profiling.sh`).

- [x] Kraken classification rates: RA_Gut has none (kreports deleted). Recovered for the other datasets from git (HANDOFF 2026-11-05); merge the recovered file on the HPC too.

- [x] **sourmash RefSeq counts are doubled** (found 2026-11-03, see `doc/sourmash_refseq_double_counting.txt`). R joins gathers to `SM_RefSeq_20250528_lineages.csv` on the bare accession number, which matches the GCA/GCF twins and versions, so `left_join` duplicates most genomes (+111–116 % on P19_Saliva GQ10–12). Counts, rarefaction depth and count-based alpha are affected. Fixed on the `merged_tables` branch (RefSeq lineages from the sylph-tax table). The relative abundances shift too (max Shannon difference 0.1–0.7 per dataset), so the sourmash RefSeq results need regenerating.

## Methods notes (for the manuscript)

- [ ] **The sample sets differ between the diversity comparisons and the statistical tests.**
  - The diversity comparisons use every sample of a dataset: alpha diversity, PCoA, pairwise distances and Procrustes.
  - The Wilcoxon tests and PERMANOVA use only the samples that have a value for the grouping variable (`grouping_variable` in `0_config.R`).
  - Maple: 88 samples in the diversity comparisons, but only the 66 with a Provenance (33 M, 33 S) in the tests.
  - The `n` in the figure labels (`dataset_names_n_lookup_alphadiv.RDS`) comes from the Wilcoxon tests, so for Maple it's the tested count.
  - State this in the methods.
- [ ] **PD mOTUs 3: 70 profiles salvaged (tell reviewers).** 70 of the 725 PD mOTUs 3 profiles (rows 1–25 plus 45 scattered rows) were empty files from a batch that died on 2025-03-29, so the merged table had 655 samples. Until 2026-11-07 the classification-rate script also gave each of them the previous sample's rate. Rerun on 2026-11-07 (job 407787, same mOTUs 3.1.0 and database). Mention the rerun in the methods or the response to reviewers, then rerun `merge_tables.sh PD` and the PD mOTUs 3 figures.
- Out of scope for this project: rarefaction curves (`9.1`/`9.2`) and the DAA scripts (`5.0`–`5.4`).
- [ ] **New datasets (Nov 2026): describe them in the methods.**
  - **Shared preprocessing:** raw paired reads were retrieved from the NCBI Sequence Read Archive. Each run's read count was checked against the archive's.
    - QC used KneadData 0.12.0, with Trimmomatic (`SLIDINGWINDOW:4:20 MINLEN:50`), then removal of host reads by mapping to the host genome with Bowtie2 (`--very-sensitive-local`).
    - These settings are the same as for Maple.
  - **Oilfield (PRJNA1196024):** oil-contaminated coastal wetland soils, surface layer (0–20 cm).
    - Paper: Commun. Earth Environ. 2025, doi:10.1038/s43247-025-02062-1.
    - 87 samples, one per site: Shengli Oilfield 44, Dagang Oilfield 43. Illumina NovaSeq 6000, 2×150.
    - **Grouping variable:** oilfield, taken from the sample names (`-SL` / `-DG`).
    - **Host removal:** human GRCh38 (contamination filter only; no host).
  - **Grapevine (PRJNA1153634):** root rot of wine grapevines, Ningxia, China.
    - Paper: Genome Biol. 2026, doi:10.1186/s13059-025-03905-y.
    - 78 samples from 39 vines: healthy vs diseased × 2 sites × 2 cultivars (Cabernet Sauvignon, Syrah) × rhizosphere soil and root. DNBSEQ-T7, 2×150.
    - Metadata decoded from the sample names (e.g. `UHDZR3`).
    - **Grouping variable:** health status or compartment (to decide).
    - **Pairing:** each vine contributes a soil and a root sample, so the compartment comparison is paired. Our tests are unpaired, so state it as a limitation, or test one compartment.
    - **Host removal:** *Vitis vinifera* PN40024 T2T reference (GCF_030704535.1).
  - **Snake (PRJNA1353898):** cloacal swabs of Eastern Indigo Snakes (*Drymarchon couperi*), Florida, 2025.
    - Paper: PLOS ONE 2026, doi:10.1371/journal.pone.0350824.
    - 80 snakes: 60 captive (one facility) and 20 wild. Illumina NovaSeq X.
    - **Grouping variable:** wild vs captive, from the BioSample location.
    - The paper splits the captive snakes by *Cryptosporidium serpentis* status (30/30), but that status isn't public.
    - **Host removal:** *Drymarchon corais* (GCA_043091225.1), the reference the authors used.
    - Status: only 3 samples are downloaded so far, as a test.

## Adding a dataset

`scripts/0_run_all.R` builds the phyloseq objects and the diversity metrics one dataset at a time. Results are cached per dataset in `Out/_Rdata/cache/<DS>/`, and the script then assembles the usual `Out/_Rdata/*.RDS` files. Only datasets without a complete cache are computed; list datasets in `rerun` to recompute them anyway.

1. HPC:
   - run the profilers;
   - add a line to `scripts/src/datasets.tsv`;
   - run `bash scripts/src/merge_tables.sh <DS>`, plus `scripts/utils/sylph_merge.sh <DS> <db> <sample tsv>` for each sylph database.
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

## To launch

- [x] **Horse on Fir** (PRJNA1102860): profiled and merged 2026-11-07, outputs in `data/Horse` on Fir (2026-11-09). Remaining: the R steps of "Adding a dataset". Original note: Merge the `merged_tables` and ILL `kneaddata_gzip` PRs first. Set `SRC` in `setup_fir.sh`, and `HORSE_FASTA` in `horse.sh`, or `HORSE_INDEX` in `setup_fir.sh` if the index is built on iv12. Steps are in `doc/HANDOFF_FIR.md`.
- [ ] **Resource benchmark** for the reviewer's runtime and memory request (written 2026-11-06, not submitted: the queue is full; scripts are local and uncommitted, held out of the `merged_tables` PR). 55 samples (5 per dataset, no Olive), all tools and databases, all on `c-3t`. Expect a few days on that node. Output: `Out/benchmark/`.
  - **Launch:** paste the blocks of `scripts/benchmark/benchmark.sh` in order, from a login or interactive shell. Everything goes through `sbatch`; no `salloc` needed. The sample block rewrites `samples.tsv` and `tools.tsv`, already built once.
  - **Pilot first:** after the sample block, set `A=1` and `K=1` to run one sample per tool. Check memory and the Kraken step below, delete `$B/jobs.tsv`, then set them back to `A=1-$(wc -l < $B/tools.tsv)` and `K=$(wc -l < $B/tools.tsv)`.
  - **The commands, once `submit` is defined:**
    ```bash
    submit bases $BM/benchmark_bases.sh $B
    submit sylph_sketch    --array=$A%10 -n 2 --mem=16G $BM/benchmark_sketch.sh sylph $B/tools.tsv $B
    SK_SYLPH=$LAST
    submit sourmash_sketch --array=$A%10 -n 1 --mem=16G $BM/benchmark_sketch.sh sourmash $B/tools.tsv $B
    SK_SM=$LAST
    # + the 3 $metaphlan generator lines, then:
    submit MPA22 --array=$A%2 $B/MPA_db2022/metaphlan.slurm.sh
    submit MPA23 --array=$A%2 $B/MPA_db2023/metaphlan.slurm.sh
    submit MPA25 --array=$A%2 $B/MPA_db2025/metaphlan.slurm.sh
    submit MOTU3 --array=$A%3 -n 6  --mem=60G  $SRC/motus3_SLURM.sh $B $B/tools.tsv
    submit MOTU4 --array=$A%2 -n 12 --mem=100G $SRC/motus4_SLURM.sh $B $B/tools.tsv
    DS=../Out/benchmark
    submit SYref --array=$A%5 -n 2 --mem=64G --dependency=afterany:$SK_SYLPH $SRC/sylph_SLURM.sh \
    	$DS $B/tools.tsv $DB_SYLPH/RefSeq_20250528_c200.syldb $DB_SYLPH/RefSeq_20250528_taxonomy.tsv.gz RefSeq
    submit SY220 --array=$A%5 -n 2 --mem=64G --dependency=afterany:$SK_SYLPH $SRC/sylph_SLURM.sh \
    	$DS $B/tools.tsv $DB_SYLPH/gtdb-r220-c200-dbv1.syldb $DB_SYLPH/gtdb_r220_metadata.tsv.gz gtdb_rs220
    submit SMref --array=$A%5 -n 1  --mem=64G --dependency=afterany:$SK_SM $SRC/gather_SLURM_fast.sh \
    	$DS $B/tools.tsv $DB_SM/RefSeq_20250528.k31.sbt.zip
    submit SM220 --array=$A%2 -n 24 --mem=64G --dependency=afterany:$SK_SM $SRC/gather_SLURM_fast.sh \
    	$DS $B/tools.tsv $DB_SM/gtdb-rs220-reps-k31.zip
    submit KB --time=48:00:00 -N 1 -c 24 --mem=200G $BM/benchmark_kraken.sh $DB_K2/k2_standard_20241228 KB "" $B/tools.tsv $K $B
    submit KG --time=48:00:00 -N 1 -c 24 --mem=700G $BM/benchmark_kraken.sh $DB_K2/k2_gtdb_genome_reps_20241109 KG _GTDB $B/tools.tsv $K $B
    ```
  - **When finished:** `check_jobs` prints nothing, then the last block runs `benchmark_gather.sh` and `scripts/10_benchmark_plot.R`, which write `resource_use.tsv`, `summary.tsv`, `benchmark_ridges.pdf` (raw CPU-hours, CPU-hours per Gbp, peak memory) and `benchmark_scaling.pdf` (CPU-hours against Gbp, per tool).
  - **Check at the pilot:**
    - Does a Kraken sample step's `MaxRSS` include the `/dev/shm` index? The R script adds the index size when the step is below it (85 and 502 GB).
    - Do the `srun` steps see the index copied by the `KB_copy` or `KG_copy` step?
    - `kraken_local.sh` runs under `nice -n10` between 06:00 and 20:00. This can slow elapsed time, not `TotalCPU`.
