# Handoff

Read at the start of a session and keep it current. Open items: `doc/TODO.md`. Older logs: `git log -p -- doc/HANDOFF.md`.

## Current state (2026-11-19, end of a Mac session)

- **Results:** `Out/_Rdata/*.RDS` from the 2026-11-12 rerun (all 12 datasets) and the 2026-11-16 beta-sensitivity objects, synced to the Mac.
- **Resource benchmark:** 55 samples (11 datasets, no Horse), all 17 tools done on cp0359 (c-3t), job ids in `Out/benchmark/jobs.tsv`. Two gaps, fixed by the HPC task below:
  - Horse is missing (the samples were picked on 2026-11-06, when Horse was on Fir only);
  - sourmash isn't comparable between databases: RefSeq is a Sequence Bloom Tree (SBT) index, which branchwater's `fastgather` can't read in any version, so it ran on 1 thread with the Python `gather`; GTDB ran on 24 threads with `fastgather`. Decided 2026-11-19: upgrade to sourmash 4.9.4 + branchwater 0.9.15, convert RefSeq to a zip, rerun the sourmash part of the benchmark only (sketch + both gathers, 24 threads each). The main analysis keeps its 4.8.11 profiles, provided the comparison in step 6 shows the same composition.
- **Mac, for the user to commit and push before the HPC session:**
  - `scripts/src/gather_SLURM_fast.sh`, `scripts/benchmark/src/benchmark_sketch.sh`: the container comes from `$SOURMASH_SIF` when set (default unchanged, 4.8.11).
  - `scripts/benchmark/benchmark.sh`: new sourmash 4.9.4 block; the old sourmash gather lines are commented out.
  - New: `scripts/benchmark/sourmash_4.9.4.def` (container recipe).
  - Stage: `git add -u && git add scripts/benchmark/sourmash_4.9.4.def && git add -f doc/HANDOFF.md doc/TODO.md`.

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
**Route B, recipe (if no plugin image):** `scripts/benchmark/sourmash_4.9.4.def` starts from `quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0` (tag confirmed) and pip-installs the plugin:
```bash
apptainer build --fakeroot $C/sourmash.4.9.4.sif scripts/benchmark/sourmash_4.9.4.def
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
  - 60 gather files in each of `$B/sm49/SM_RefSeq_20250528/` and `$B/sm49/SM_gtdb-rs220-reps/`;
  - every row of `sourmash_sketch`, `SMref` and `SM220` has a COMPLETED job from the new submission (`10_benchmark_plot.R` keeps the last COMPLETED job per label × row, so a failed new task would silently leave the 4.8.11 row);
  - every other tool has rows 56–60.
- Record the node's hardware for the methods: `scontrol show node cp0359 | grep -E "CPUTot|RealMemory"; ssh cp0359 lscpu | grep "Model name"`. Put it in this file.

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
