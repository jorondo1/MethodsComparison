# Handoff

Read at the start of a session and keep it current. Open items: `doc/TODO.md`. Older logs: `git log -p -- doc/HANDOFF.md`.

## Current state (2026-11-19, end of a Mac session)

- **Results:** `Out/_Rdata/*.RDS` from the 2026-11-12 rerun (all 12 datasets) and the 2026-11-16 beta-sensitivity objects, synced to the Mac.
- **Resource benchmark:** 60 samples since 2026-11-18 (Horse added on the HPC; see step 3). Before that, 55 samples (11 datasets, no Horse), all 17 tools done on cp0359 (c-3t), job ids in `Out/benchmark/jobs.tsv`. Two gaps, fixed by the HPC task below:
  - Horse is missing (the samples were picked on 2026-11-06, when Horse was on Fir only);
  - sourmash isn't comparable between databases: RefSeq is a Sequence Bloom Tree (SBT) index, which branchwater's `fastgather` can't read in any version, so it ran on 1 thread with the Python `gather`; GTDB ran on 24 threads with `fastgather`. Decided 2026-11-19: upgrade to sourmash 4.9.4 + branchwater 0.9.15, convert RefSeq to a zip, rerun the sourmash part of the benchmark only (sketch + both gathers, 24 threads each). The main analysis keeps its 4.8.11 profiles, provided the comparison in step 6 shows the same composition.
- **Mac, for the user to commit and push before the HPC session:**
  - `scripts/src/gather_SLURM_fast.sh`, `scripts/benchmark/src/benchmark_sketch.sh`: the container comes from `$SOURMASH_SIF` when set (default unchanged, 4.8.11).
  - `scripts/benchmark/benchmark.sh`: new sourmash 4.9.4 block; the old sourmash gather lines are commented out.
  - The container recipe is `$ILL_PIPELINES/containers/sourmash.4.9.4.def` (moved there on 2026-11-19, next to the other sourmash recipes and the .sif; not in this repo).
  - Stage: `git add -u && git add -f doc/HANDOFF.md doc/TODO.md`.

## Task for the HPC session: Horse + sourmash 4.9.4 in the resource benchmark

`git pull` first. Steps 1–2 build the tools, 3 adds Horse, 4 submits, 5–6 check. Keep every benchmark job on c-3t (`submit` does it).

### 1. Build the sourmash 4.9.4 container (and learn how)

**How Apptainer containers work.** A `.sif` file is a read-only, compressed filesystem image: you can't upgrade the software inside it. To "update" a container, you build a new image and give it a new name (`sourmash.4.9.4.sif`), and keep the old one, because the published profiles were made with it. Two ways to build:
- **Pull** an existing Docker image (`apptainer build new.sif docker://...`): Apptainer downloads its layers and squashes them into one SIF. No root needed. BioContainers (quay.io/biocontainers) publishes an image for every Bioconda package.
- **Build from a recipe** (a `.def` file): start from an image (`Bootstrap: docker`, `From: ...`), then run install commands in `%post`, check in `%test`, and document in `%labels`. Running `%post` needs root inside the build, which `--fakeroot` provides if the cluster allows it.

Useful commands:
- `apptainer inspect --deffile old.sif` prints the recipe an image was built from (empty for a pulled image);
- `apptainer inspect --labels old.sif` shows its labels;
- `apptainer exec old.sif <cmd>` runs a command inside it.

**Do it on a compute node** (building is CPU- and disk-heavy): `salloc -c 2 --mem=8G --time=2:00:00 -A def-ilafores`. Then:
```bash
module load StdEnv/2020 apptainer/1.1.5     # or a newer apptainer: module spider apptainer
export APPTAINER_CACHEDIR=$SLURM_TMPDIR/apptainer_cache APPTAINER_TMPDIR=$SLURM_TMPDIR  # keep the layer cache out of $HOME (quota)
C=$ILAFORES/programs/ILL_pipelines/containers   # where the other containers live; any folder works if SOURMASH_SIF in benchmark.sh follows
# a) How was the current container built, and does it have branchwater?
apptainer inspect --deffile $C/sourmash.4.8.11.sif
apptainer exec $C/sourmash.4.8.11.sif sourmash info -v      # lists plugins, branchwater among them
```
**Route A, pull (try first):** look for a BioContainers image of the plugin. It depends on sourmash, so it brings the sourmash CLI along:
```bash
curl -s "https://quay.io/api/v1/repository/biocontainers/sourmash_plugin_branchwater/tag/?limit=100&onlyActiveTags=true" \
	| grep -o '"name": *"[^"]*"' | sort -u        # look for 0.9.15--<build>
apptainer build $C/sourmash.4.9.4.sif docker://quay.io/biocontainers/sourmash_plugin_branchwater:<that tag>
```
**Route B, recipe (used: there is no BioContainers image of the plugin).** `$ILL_PIPELINES/containers/sourmash.4.9.4.def` follows `sourmash.4.8.11.def` (Ubuntu, apt, pip compiles the plugin), with Ubuntu 24.04 (the plugin needs Python ≥ 3.11) and Rust 1.85 from rustup (the plugin needs Rust ≥ 1.85). PyPI has only the plugin's source, so it compiles: about 270 Rust crates.
```bash
cd $C && apptainer build --force --fakeroot sourmash.4.9.4.sif sourmash.4.9.4.def
```
If `--fakeroot` is refused, stop and ask the user; don't improvise a workaround.

**Check the new image** (either route):
```bash
S=$C/sourmash.4.9.4.sif
apptainer exec $S sourmash --version          # 4.9.4. If a pull gave another 4.9.x, name the file after it and update SOURMASH_SIF in benchmark.sh
apptainer exec $S sourmash info -v            # branchwater 0.9.15
apptainer exec $S sourmash scripts fastgather --help   # note the defaults of --threshold-bp and --cores, for step 6
```

### 2. Convert the RefSeq SBT to a zip collection

**Done 2026-11-19 (job 408976, `sig cat`, 4.8.11 container):** `$ILAFORES/ref_dbs/sourmash_db/RefSeq_20250528.k31.zip`, 1.9 GB. `sig fileinfo` matches the SBT: 67,584 sketches, DNA k=31 scaled=1000, 226,849,503 hashes. The "1,649 MinHashes encountered multiple times" are identical sketches (probably GCA/GCF twins), all kept. The commands are logged in `scripts/src/sourmash_refseq_db.sh` (Index section). Index sizes of every database: `Out/benchmark/database_sizes.tsv`.

A zip is branchwater's recommended format. `sig cat` should read the SBT's leaves (all sourmash commands load SBTs, but the docs don't say so for `sig cat` specifically), so check the counts:
```bash
DB_SM=$ILAFORES/ref_dbs/sourmash_db; S=$ILAFORES/programs/ILL_pipelines/containers/sourmash.4.9.4.sif
sbatch -A def-ilafores -c 1 --mem=64G --time=12:00:00 -J sm_refseq_zip -o logs/sm_refseq_zip-%j.out --wrap \
	"module load StdEnv/2020 apptainer/1.1.5; singularity exec -B $ILAFORES $S sourmash sig cat $DB_SM/RefSeq_20250528.k31.sbt.zip -o $DB_SM/RefSeq_20250528.k31.zip"
# then compare: same number of sketches, ksize 31, scaled 1000, abundances, and total hashes
singularity exec -B $ILAFORES $S sourmash sig fileinfo $DB_SM/RefSeq_20250528.k31.sbt.zip
singularity exec -B $ILAFORES $S sourmash sig fileinfo $DB_SM/RefSeq_20250528.k31.zip
```
If `sig cat` can't read the SBT, build the zip from the original per-genome sketches instead. `scripts/src/sourmash_refseq_db.sh` wrote them to `/fast2/def-ilafores/refseq_genomes/20250528_signatures/`. A text file listing their paths, one per line, is itself a valid input: `find ... -name '*.sig' > list.txt; sourmash sig cat list.txt -o ...zip` (add `-B /fast2/def-ilafores`). Keep the `.sbt.zip`. The lineage file stays the same.

### 3. Append 5 Horse samples (rows 56–60)

**Done on the HPC 2026-11-18** (before this plan reached the HPC), so skip this step. How it differed from the plan:
- Horse's post-QC reads were on Fir only, and Fir was down. The 5 runs were downloaded from the SRA (`data/Horse/raw_hpc/`) and kneaddata was rerun on the HPC, into `data/Horse/preproc_hpc/`. All 18 READ COUNT lines of each log are identical to Fir's.
  - Built on the HPC for this: the EquCab3 index `$FAST/host_genomes/EquCab3_index/EquCab3`, and the adapter copy `$FAST/adapters/TruSeq3-PE-2.fa` (see the TODO adapter item).
  - cpubase nodes cp04xx can't run StdEnv/2023 (x86-64-v3) builds ("Illegal instruction"). Use `-p bio` for those jobs.
- Picked with the benchmark.sh rule, ranked by kneaddata's final pair1 count, not the paired_1 size (no reads to stat then): SRR28771685, SRR28795704, SRR28772227, SRR28772014, SRR28771816 (0.5–13M pairs). Backups: `*.bak_pre_horse` (rather than `.bak_55`).

**Plan text, kept for reference:**

**Don't rerun the sample block of `benchmark.sh`:** it rewrites `samples.tsv` and `tools.tsv`, which would undo the 2026-11-12 nfs-ip34 path fix in `tools.tsv`.
1. Check that Horse's post-QC reads are on the HPC, under `$ILAFORES` (the Kraken container sees only `$ILAFORES`): `head -2 $MC/data/Horse/preproc/preprocessed_reads.sample.tsv`, then `ls` those paths. On Fir only → stop and ask.
2. Paste the preamble of `benchmark.sh` (variables, `N`, `MIN_MB`) and the `submit()` function, then:
```bash
cp $B/samples.tsv $B/samples.tsv.bak_55; cp $B/tools.tsv $B/tools.tsv.bak_55; cp $B/bases.tsv $B/bases.tsv.bak_55
grep -P '^Horse\t' scripts/src/datasets.tsv | while IFS=$'\t' read -r ds tsv; do
	tsv=$(eval echo "$tsv")
	while IFS=$'\t' read -r sample p1 p2 u1 u2 _; do
		[[ -s $p1 && -s $p2 && -s $u1 && -s $u2 ]] && (( $(stat -c %s "$p1") >= MIN_MB * 1024 * 1024 )) &&
			printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(stat -c %s "$p1")" "$sample" "$p1" "$p2" "$u1" "$u2"
	done < "$tsv" | sort -n | awk -F'\t' -v OFS='\t' -v n=$N -v ds=$ds '
		{ row[NR] = $0 }
		END { for (i = 1; i <= n && i <= NR; i++) { k = (NR <= n) ? i : int(1 + (i - 1) * (NR - 1) / (n - 1) + 0.5)
			split(row[k], f, "\t"); print ds, f[2], f[3], f[4], f[5], f[6], f[1] } }'
done > $B/samples_horse.tsv
cat $B/samples_horse.tsv                      # 5 rows, paths under $ILAFORES
cat $B/samples_horse.tsv >> $B/samples.tsv
cut -f2-6 $B/samples_horse.tsv >> $B/tools.tsv
wc -l $B/samples.tsv $B/tools.tsv             # 60 each
submit bases_horse -n 4 --mem=8G --time=12:00:00 -A def-ilafores \
	--wrap "ml seqkit; cut -f3-6 $B/samples_horse.tsv | tr '\t' '\n' | xargs seqkit stats -T -j 4 | tail -n +2 >> $B/bases.tsv"
```

### 4. Submit

**sourmash 4.9.4 pilot, rows 1 and 30 (2026-11-19):** sourmash_sketch 408977, SMref 408978, SM220 408979 (24 cores, 64G; fastgather defaults checked: `-t` 50 kb like gather, `-c` all available, so the 24 comes from `-n 24`). All COMPLETED. Peak memory 2.5G (64G is ample); row 30 (PD, 52 Gbp): sketch 2 h 14, SMref 21 min wall / 6.2 CPU-h, SM220 22 min / 5.1 CPU-h. Why 24 cores: `CPU_h` is sacct TotalCPU (CPU time used), so it barely depends on the thread count; 24 matches Kraken and 4.8.11 GTDB.

**sourmash 4.9.4 full run, rows `2-29,31-60` (2026-11-19):** sourmash_sketch 408985, SMref 408986, SM220 408987 (same block, `A=2-29,31-60`; the pilot rows were left out because the gather script skips rows whose output exists, so a rerun would record a near-zero time). When they finish: step 5 below, then the version check.
- **All done 2026-11-20:** sketch, SMref, SM220 58/58 COMPLETED; 56 RefSeq + 60 GTDB gather files; `sacct.tsv` regenerated (3,344 records). Last COMPLETED job per label × row: every label has 60 rows, and every sourmash row points to a 4.9.4 job (408977–79 for rows 1 and 30, 408985–87 for the rest).
- **Node hardware (cp0359, for the methods):** 4 × Intel Xeon E7-8870 v3 @ 2.10 GHz (Haswell), 72 cores (18 per socket, no hyper-threading), 2.83 TiB RAM (RealMemory 2,896,361 MB).
- **Memory, RefSeq > GTDB with 4.9.4 too** (median peak 0.85 G against 0.51 G; all < 1.5 G): fastgather keeps every prefetch candidate (≥ 50 kbp overlap) in memory, and RefSeq's within-species redundancy (strains, GCA/GCF twins) gives 2–4× more candidates on human samples (mean 5,970 against 2,913). Where candidate counts are close (row 50: 18.4k against 12.4k), the memory is the same (1.47 G against 1.50 G).
- **SMref Horse rows 56, 57, 59, 60 have no gather file, and that's expected:** no RefSeq genome reaches the 50 kbp overlap, so fastgather stops after prefetch ("No search signatures loaded, exiting") and writes nothing. 4.8.11 did the same on these 4 rows ("Prefetch found 0 signatures with overlap >= 50.0 kbp"). Row 58 has 1 match in both versions. The search ran in full, so the timings stand. Expect **56** RefSeq gather files, not 60. The `mv` error in those logs comes from the missing tmp file. The tasks still show COMPLETED because `gather_SLURM_fast.sh` doesn't propagate the exit code. `compare_sourmash_versions.R` keeps only samples with a 4.9.4 file (line 41), so these 4 drop out of the RefSeq comparison (empty in both versions).

**sourmash version check (task 2):** `bash scripts/benchmark/sourmash_version_check.sh` (sketch md5s, RefSeq lineages, tax summaries of the 4 gather sets with the 4.8.11 tax code → `$B/version_check/`), then `module load StdEnv/2023 r/4.4.0; Rscript scripts/benchmark/compare_sourmash_versions.R` (→ `$B/sourmash_version_check.tsv` + verdict per database). Pilot result (rows 1, 30):
- Sketches: identical md5.
- GTDB: **identical** (same genomes, order and `f_unique_weighted`).
- RefSeq: **equivalent**, not identical. Of 198–342 genomes per sample, ~9 differ. They are ties: genomes with the same overlap and `f_unique_weighted`, broken differently by the SBT (4.8.11 Python gather) and the zip (Rust fastgather). Mostly within one species, so no effect. Where the tied genomes carry different species names (PD: *Faecalibacterium taiwanense* ↔ *F.* sp. i25-0019-C1, 1.1% of the hashes), the abundance moves between them. Bray-Curtis 0.003 and 0.025; Jaccard ≥ 0.99; top 10 same.
- Anchor: the 4.8.11 benchmark summary equals the production merged table (RefSeq, 2/2).

**Full result (2026-11-20, all 60 rows):**
- Sketches: 60/60 identical.
- **GTDB: identical** in all 60 samples (same genomes, order, `f_unique_weighted`; Bray-Curtis 0). Anchor 20/20.
- **RefSeq: equivalent.** The script prints "DIFFERENT" because its rule needs the top 10 to match exactly, but every difference is a tie:
  - 56 samples (the 4 Horse rows without matches drop out). Bray-Curtis median 0.0009, max 0.043; Jaccard min 0.81; genome lists identical in 5 samples only. As in the pilot, the 4.8.11 SBT and the 4.9.4 zip break ties between genomes with equal overlap differently.
  - Top 10 differs in 3 samples, each a swap at the rank 10/11 boundary between species whose abundances differ by ≤ 0.0002 (GQ13: *S. australis* in, *N.* sp. RH3002v2g out; SRR5091560: *B. fragilis* ↔ *B. cellulosilyticus*, 0.0172 against 0.0171; Moss S-18-POLPIL-B: *Methylobacterium* sp. NMS14P ↔ *Bradyrhizobium lablabi*, all < 0.0002 of the hashes).
  - Highest Bray-Curtis, Moss S-18-POLPIL-B (0.043): a sample with barely any RefSeq signal (top species 0.16 % of the hashes), so a few tied genomes weigh a lot after renormalising.
  - Largest count change (476k, PD SRR19064482): the *Faecalibacterium taiwanense* ↔ *F.* sp. i25-0019-C1 tie seen in the pilot.
  - **Conclusion: the main analysis keeps its 4.8.11 profiles.**
- Anchor RefSeq 54/56. The 2 misses (GQ19, GQ38) don't come from the upgrade: their production sketches are from older read files (2024 `.fastq` on `/nfs3_ib/...`) than the current `preproc/*.fastq.gz` the benchmark sketched (GQ19: 205,897 against 209,816 hashes; GQ38 4 % larger in production). The other 6 P19 sketches have identical sizes. See TODO.

**Horse rows, every tool: done 2026-11-18** (jobs 408887–408903, all COMPLETED; the 30 Horse kreports have reads). Every tool has 60 rows in `resource_use.tsv`. sourmash 4.8.11 also ran on rows 56–60 (sourmash_sketch 408889, SMref 408900, SM220 408901). The 4.9.4 rerun of all 60 rows below supersedes them, since the plot keeps the last COMPLETED job. The sourmash 4.9.4 part is submitted (pilot + full run above).

- **Horse rows, every tool except sourmash:** set `A=56-60` and paste the `benchmark.sh` blocks: the sylph sketch (**skip the `sourmash_sketch` line of that block**), MetaPhlAn, mOTUs, SYref/SY220. Then Kraken, with the row list as 7th argument:
  ```bash
  submit KB --time=48:00:00 -N 1 -c 24 --mem=200G $BM/benchmark_kraken.sh $DB_K2/k2_standard_20241228 KB "" $B/tools.tsv 60 $B 56,57,58,59,60
  submit KG --time=48:00:00 -N 1 -c 24 --mem=700G $BM/benchmark_kraken.sh $DB_K2/k2_gtdb_genome_reps_20241109 KG _GTDB $B/tools.tsv 60 $B 56,57,58,59,60
  ```
  - MetaPhlAn: don't rerun the generator lines; submit the existing `$B/MPA_db20*/metaphlan.slurm.sh`. First check that they read their row from `$B/tools.tsv` by `SLURM_ARRAY_TASK_ID`.
  - mOTUs 4 needed 400G on the largest PD sample (52 Gbp): give a deep Horse sample more than 100G.
- **sourmash 4.9.4, all 60 rows:** set `A=1-60` and paste the new block "sourmash 4.9.4 + branchwater" of `benchmark.sh`. It exports `SOURMASH_SIF`, writes to `$B/sm49/` (new sketches too; the 4.8.11 outputs stay in `$B/SM_*` and `$B/signatures`), and reuses the labels `sourmash_sketch`, `SMref`, `SM220`.

### 5. When everything has finished

- `bash $BM/benchmark_gather.sh $B` (regenerates `sacct.tsv`).
- Check:
  - 56 gather files in `$B/sm49/SM_RefSeq_20250528/` (4 Horse rows have no RefSeq match, see step 4) and 60 in `$B/sm49/SM_gtdb-rs220-reps/`;
  - every row of `sourmash_sketch`, `SMref` and `SM220` has a COMPLETED job from the new submission (`10_benchmark_plot.R` keeps the last COMPLETED job per label × row, so a failed new task would silently leave the 4.8.11 row);
  - every other tool has rows 56–60.
- Record the node's hardware for the methods: done, see step 4 (ssh to cp0359 is refused; an `sbatch -w cp0359 --wrap lscpu` works, with `-o` on shared storage, not `/tmp`).

### 6. Sub-agent: does 4.9.4 change the taxonomic composition?

Spawn a sub-agent for this. Give it this brief:
- **Question:** for the 55 samples profiled by both versions, does the species composition differ between the sourmash 4.8.11 outputs and the 4.9.4 ones, for each database?
  - RefSeq: Python `gather` on the SBT against `fastgather` on the zip.
  - GTDB r220 reps: `fastgather` 4.8.11 against 4.9.4.
- **Inputs:**
  - 4.8.11: `Out/benchmark/SM_RefSeq_20250528/*_gather.csv` and `Out/benchmark/SM_gtdb-rs220-reps/*_gather.csv`;
  - 4.9.4: the same names under `Out/benchmark/sm49/`;
  - samples: `Out/benchmark/samples.tsv`.
- **Method:**
  - Apply the production taxonomy step to both: `merge_sourmash()` in `scripts/src/merge_tables.sh`, i.e. `scripts/src/sourmash_tax_summary.py` with the 4.8.11 container and the same lineage files. Then parse to species as `scripts/utils/parse_merged.R` does: same abundance column, same RefSeq lineage handling (the double-counting fix).
  - `fastgather` renames columns (`match_name`, `match_md5`, `match_filename`; no `potential_false_negative`): check that the tax step reads both versions.
  - Compare the threshold defaults first: `gather` uses 50 kbp; for `fastgather`, see its `--help` in step 1.
  - Light work: the interactive node (≤ 2 cores) or a small `salloc`.
- **Per sample × database:**
  - species richness;
  - Jaccard of the species sets;
  - Bray-Curtis on relative abundances;
  - Shannon difference;
  - fraction of the sample classified;
  - number of genomes matched;
  - species found by only one version, with their abundance.
- **Output:**
  - `scripts/benchmark/compare_sourmash_versions.R` (tracked) and `Out/benchmark/sourmash_version_check.tsv`;
  - a short verdict per database, in this file: for example, "unchanged" if Bray-Curtis ≤ 0.05 and the top 10 species are the same in every sample; otherwise the samples that differ and why.
- **If composition changed**, ask the user before anything else: the main analysis's sourmash RefSeq (and maybe GTDB) profiles might then need a rerun with 4.9.4.

### 7. Hand back to the Mac

Update this file (hardware, verdict, anything unexpected). On the Mac: `bash scripts/utils/sync_from_hpc.sh <user@host>`, then rerun `scripts/10_benchmark_plot.R` and update `scripts/benchmark/benchmark_methods.txt`:
- Horse in, 60 samples, the new Gbp range and median;
- sourmash 4.9.4 with 24 threads for both databases (the main analysis used 4.8.11);
- the node's hardware.

## Sync to the Mac

From the Mac repo root: `bash scripts/utils/sync_from_hpc.sh <user@host>` (`--dry-run` to preview; key from `$SSH_KEY`, default `~/.ssh/ronj`). It pulls `Out/classification_rates/*.tsv`, the `Out/_Rdata/*.RDS` the figures read, and `Out/benchmark/*.{tsv,pdf}`. `-u` keeps newer local files.

## Next steps on the Mac (not blocked by the HPC)

- Finish `2_alpha_compare.R` (triad figures; it writes `Out/Manuscript/PLOS/Fig1.tiff`, so set the May version aside first).
- Render `3_`–`5_`, `6.1`, `6.2` with the 2026-11-12 objects.
- Open: how sylph appears in the figures. Out of scope: `6._FigX_TaxAssignment.R`, `9.x`, DAA `5.x`.

## Decisions on record

- **Samples:** every post-QC sample enters the diversity comparisons; samples without metadata only drop out of the Wilcoxon tests and PERMANOVA. `scripts/src/dropped_samples.tsv` removes samples from everything (P19_Gut controls and mock, PD's 5 SP/SC runs, AD_Skin's 78 runs without metadata). 2,216 samples post-QC, 2,185 with metadata.
- **Rarefaction** (2026-11-16): at the minimum sample, per dataset × database; no depth floor. A sensitivity analysis with an outlier rule (< 10 % of the table median) may come later.
- **Benchmark:** every job on cp0359 (c-3t), so the hardware is identical. Cost model: Gamma GLM, identity link (variance grows with the mean; `a` = overhead per sample). Kraken's index copy (< 0.15 CPU-h per dataset) is left out.
- **Planned:** a supplementary sample table (per sample: post-QC reads, group, why it was dropped, presence per database after rarefaction).

## HPC cleanup (the user's call)

- Stale: `Out/_Rdata/{ps_filt,ps_radEmu,ps_taxranks,ps_rare_taxranks}` (2025), `data/Grapevine/_partial_*`, `data/Mosquito/`.
- Jobs 405784_27–30 are stuck COMPLETING on cp0419; only admins can clear them.
- `Out/benchmark/`: outputs of the old plot script (`benchmark_cost.pdf`, `cost_projection.tsv`, `benchmark_fit.pdf`, `benchmark_ridges.pdf`, `benchmark_scaling.pdf`, `summary.tsv`), on both machines.

## Rules

- No commits or pushes unless asked. Never `git reset --hard` or `git clean`. Sync clones with fetch + mixed reset.
- Don't rerun profilers unless asked (the benchmark reruns above are asked). Never delete per-sample profiler outputs.
- Never replace a script that running jobs read (NFS: they stop silently with exit 0).
- HPC interactive node: at most 2 cores unless the user grants more. Heavy work goes through sbatch.
- Out of scope: GTDB 214, GenBank, full GTDB rs220, Olive.
