# Fir handoff: resource benchmark

You are on Fir and have no memory from the HPC sessions. Read this first, then `CLAUDE.md`. This is the last planned task on Fir: the resource benchmark (CPU time and peak memory per tool × database), for a reviewer's runtime and memory request. Horse is finished (see the end of this file).

## Layout (`scripts/fir/fir_env.sh`; source it in every session)
| Variable | Path | Holds |
|---|---|---|
| `REPO` | `/project/6047190/ronj2303/MethodsComparison` | git repo. Check out branch `benchmark` |
| `MC` | `~/scratch/mc_work` | working copy: jobs run with `-D $MC`; `scripts` is a symlink to `$REPO/scripts` |
| `ILL_PIPELINES` | `~/scratch/ILL_pipelines` | branch `kneaddata_gzip`; containers in `containers/` |
| `DBS` | `~/scratch/ref_dbs` | every database (`metaphlan4_db/`, `sourmash_db/`, `sylph_db/`, `kraken2_dbs/`, `motus3/`) |
| `ILAFORES` | `~/scratch` | also holds `programs/conda_envs/sylph` and `programs/motus3_env` |

Ignore `~/scratch/MethodsComparison`: it's an old clone from February.

## Benchmark design
- **60 samples:** 5 per dataset × 12 datasets. Eleven sets of 5 were picked on the HPC (`samples.tsv`); Horse's 5 are picked here, from reads already on Fir. Each dataset's samples are evenly spaced by read file size and skip those under 10 MB.
- **Every tool and database,** as in the manuscript: MetaPhlAn 2022/2023/2025, mOTUs 3 and 4, sylph (RefSeq, GTDB 220), sourmash (RefSeq, GTDB 220 reps), Kraken + Bracken (standard and GTDB, confidences 0.10, 0.45 and 0.90).
- **One Slurm record per sample × tool:** one array task per sample. Sketching (sylph, sourmash) runs as its own array, so its cost is kept apart from profiling. Kraken copies its index into /dev/shm once per index, then runs one `srun` step per sample × confidence.
- **Same hardware:** this is why the benchmark runs on Fir, where the CPU nodes are all the same type. On the HPC it would have needed the single c-3t node.
- **Outputs** in `$MC/Out/benchmark/`:
  - `samples.tsv`, `tools.tsv`, `bases.tsv` (seqkit stats per read file);
  - `jobs.tsv` (label → job id; the last id of a label wins);
  - `sacct.tsv`, from `benchmark_gather.sh`.
  - `scripts/10_benchmark_plot.R` turns these into `resource_use.tsv`, `summary.tsv`, `benchmark_ridges.pdf` and `benchmark_scaling.pdf`. Run it on the HPC or the Mac: R on Fir may lack ggridges and patchwork. Bring back the `.tsv` files only.

## Files (branch `benchmark`)
- `scripts/benchmark/benchmark.sh`: copy-paste blocks, written for the HPC.
- `scripts/benchmark/src/`:
  - `benchmark_bases.sh` (seqkit stats);
  - `benchmark_sketch.sh` (sylph or sourmash sketch, one sample);
  - `benchmark_kraken.sh` (index into /dev/shm, then one `srun` step per sample × confidence);
  - `benchmark_gather.sh` (sacct records).
- The tool scripts are the usual ones in `scripts/src/`, the same ones that ran Horse here: `motus3_SLURM.sh`, `motus4_SLURM.sh`, `sylph_SLURM.sh`, `gather_SLURM_fast.sh`, plus ILL's MetaPhlAn generator.

## Reads
The user rsyncs them from the HPC into `$ILAFORES/benchmark_reads/` (about 154 GB). The HPC side is a folder of symlinks, `Out/benchmark/reads/<dataset>/<sample>/`, and `rsync -L` copies the real files:
```bash
source /project/6047190/ronj2303/MethodsComparison/scripts/fir/fir_env.sh
rsync -avL --partial --no-g --no-p \
	ronj2303@hpc-bio.ccs.usherbrooke.ca:/net/nfs-bio/jbod2/def-ilafores/analysis/MethodsComparison/Out/benchmark/reads/ \
	$ILAFORES/benchmark_reads/
```
`samples.rel.tsv`, in the same folder, is `samples.tsv` with paths relative to `benchmark_reads/`. Columns: dataset, sample, paired_1, paired_2, unmatched_1, unmatched_2, paired_1 bytes.

## What to change in `benchmark.sh` before running it here
Work on a copy, or edit before anything is submitted. Never edit a script that running jobs read (see Lessons). The `src/` scripts are passed to `sbatch` directly, so they are spooled and safe to edit between submissions.

1. **Environment block:**
   - Source `fir_env.sh` instead of setting `MC=$ILAFORES/analysis/...`.
   - Set `DB_MPA=$DBS/metaphlan4_db`, `DB_SM=$DBS/sourmash_db`, `DB_SYLPH=$DBS/sylph_db`, `DB_K2=$DBS/kraken2_dbs`.
   - `B=$MC/Out/benchmark` must stay under `$MC`: the sylph and sourmash calls pass `DS=../Out/benchmark`, a path relative to `$MC/data/`.
2. **Sample block:** replace the picker (it reads `datasets.tsv` and HPC paths) with:
   ```bash
   R=$ILAFORES/benchmark_reads
   awk -F'\t' -v OFS='\t' -v r=$R '{ for (i = 3; i <= 6; i++) $i = r "/" $i; print }' $R/samples.rel.tsv > $B/samples.tsv
   ```
   Then append 5 Horse samples from `$MC/data/Horse/preproc/preprocessed_reads.sample.tsv`. Keep the original picker's rule: skip a sample with a missing read file or `paired_1` under `MIN_MB`, sort by `paired_1` size, take 5 evenly spaced rows. Reuse its `while … | sort -n | awk` pipeline, with `ds=Horse` and that tsv. The output has the same 7 columns.
   - Then `cut -f2-6 $B/samples.tsv > $B/tools.tsv`, as in the script.
   - Check: 60 rows, unique sample IDs, and every file exists.
3. **`submit()`:**
   - Drop `-p c-3t`. Keep `-D $MC -o -e`: the tool scripts hardcode HPC `#SBATCH -D/-o`, so these overrides are required.
   - Fir needs a `--time` on every call. Horse's limits worked: mOTUs 6 h, MetaPhlAn 12 h, sylph 3 h, sourmash 12 h, sketch 12 h, Kraken 48 h.
   - The `%N` array throttles were set for a busy shared node. They can go up here. Jobs still share nodes, which adds a little noise to elapsed time but barely affects `TotalCPU`.
4. **`benchmark_sketch.sh`:** the sourmash container path is hardcoded as `$ILAFORES/programs/ILL_pipelines/containers/...`. On Fir it's `$ILL_PIPELINES/containers/sourmash.4.8.11.sif`. Change it to `${ILL_PIPELINES:-$ILAFORES/programs/ILL_pipelines}/containers/...`. The sylph env path (`$ILAFORES/programs/conda_envs/sylph`) already matches.
5. **`benchmark_bases.sh`:** `ml seqkit` may need `module load StdEnv/2023 seqkit`. Check with `module spider seqkit`.
6. **Kraken (`benchmark_kraken.sh`):**
   - **Memory:** `--mem` covers the index in /dev/shm: standard about 85 GB, GTDB about 502 GB. Horse's whole-node GTDB job reached 633 GB MaxRSS here. Keep `--mem=700G` for KG, or use `--mem=0` on a whole node.
   - **If /dev/shm is too small,** fall back to `$SLURM_TMPDIR` as `scripts/fir/kraken_node.sh` does. The plot script adds the index size to a sample step's memory when the step reports less, with the sizes in `index_gb`. If the index is not in /dev/shm, set `index_gb` to 0 in `10_benchmark_plot.R`.
   - **Threads:** keep `-c 24`, the per-step thread count, so Kraken runs with the same threads as on the HPC plan. One job per index, all samples one after the other: expect about a day for KG.
7. **MetaPhlAn generator:** same call as in `scripts/fir/horse.sh`. `--slurm_log $B/logs` and the walltime are already arguments.
8. **`watch.sh`:** it only watches job names matching `NAMES`. The benchmark jobs are named `bm_<label>`, so add `bm_` to `NAMES`, or pass a new pattern.

## Order of work
1. In `$REPO`: `git fetch && git switch benchmark`. The user commits and pushes; don't commit yourself.
2. Wait for the reads rsync. Check: 60 × 4 files, none empty.
3. Make the changes above, then paste the blocks of `benchmark.sh` one at a time.
4. **Pilot:** `A=1`, `K=1`. One sample per tool. Check that each job wrote its output, and look at `sacct -j <id> --format=JobID,JobName,State,Elapsed,TotalCPU,MaxRSS -P`. For Kraken, check that the sample steps see the index copied by the `KB_copy`/`KG_copy` step.
   - Then delete `$B/jobs.tsv` and the pilot outputs, set `A=1-$(wc -l < $B/tools.tsv)` and `K=$(wc -l < $B/tools.tsv)`, and submit everything.
5. When `check_jobs` prints nothing:
   - `bash $BM/benchmark_gather.sh $B`.
   - Check `sacct.tsv` for FAILED, OUT_OF_MEMORY or TIMEOUT rows, and rerun those labels. The last job id of a label is the one plotted.
   - Copy `samples.tsv`, `bases.tsv`, `jobs.tsv` and `sacct.tsv` back to `Out/benchmark/` on the HPC, and plot there.
6. Update this file and `doc/HANDOFF.md` (tag [Fir]).

## Rules from the user
- Don't commit or push unless asked. Never `git reset --hard` or `git clean`.
- Don't rerun profilers unless asked. Never delete per-sample profiler outputs.
- Light commands only on login nodes; heavy work goes through `sbatch`.
- Don't spend time on the Olive dataset, GTDB 214 or GenBank.
- Report outcomes plainly, including anything skipped or untested.

## Lessons
- **Never replace a script that running jobs may be reading,** even with an atomic rename. On NFS, bash on the compute node can stop silently, with exit 0: the job shows `COMPLETED` and has no output. Scripts passed directly to `sbatch` are spooled and safe; scripts they call are not. `COMPLETED` alone is not proof: compare outputs against tasks.
- **`gather_SLURM_fast.sh`** writes `*_gather.tmp.csv` and renames it on success. Once, an empty `SIG` made `sketch` overwrite a read file (fixed). In the pilot, check that the reads are intact.
- **`/dev/shm`** counts against job memory and is private to a job.
- **`sbatch --export=ALL,VAR=x`:** such jobs were cancelled by root on the HPC. Pass script arguments.
- **Kraken:** a kreport left by a killed run is reused as if finished, so Bracken runs on a truncated report. After any kill, move aside the kreports that have no `_bracken/*_bracken_S.MPA.TXT`.
- **mOTUs 3** needs `python/3.8.10` and `/cvmfs/soft.mugqic` (bwa). It worked for Horse.

## Horse (done)
- 744 samples, kneaddata with human + horse (EquCab3), all profilers. Merged on 2026-11-07 and rsynced to `$REPO/data/Horse`, then to the HPC.
- Sparse because RefSeq covers horse gut poorly: sylph RefSeq 147 samples, KB90 734, sourmash RefSeq 205.
- Reads: `$MC/data/Horse/preproc` (397 GB) and raw reads in `~/scratch/Horse/raw`. Keep `preproc` until the benchmark's Horse samples have run. After that, deleting them is the user's call.
