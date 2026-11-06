# Horse on Fir: handoff for the Claude session there

You are on Fir (or another Alliance cluster) and have no memory from the HPC sessions. Read this first, then `CLAUDE.md`, `doc/TODO.md` and `doc/HANDOFF.md`. Your first job is to monitor the processing of the Horse dataset for the user.

## Goal
Run Horse (PRJNA1102860, feces metagenome, shallow shotgun) through kneaddata and all profilers, as an extra dataset for the methods comparison. The user supplies an accession list and picks 500–700 of the 2,778 runs. Median run: 1.05 Gbp. The project also has 30 `metagenome` and 30 `synthetic metagenome` runs, likely controls.

## Layout
| Path | Holds |
|---|---|
| `~/scratch` (`$ILAFORES`) | `programs/`, `ref_dbs/`, `fast/`, `MethodsComparison/` (`$MC`). All jobs run here |
| `/project/def-ilafores/ronj2303/Horse/raw` (`$RAW`) | Raw reads. Same place as `~/projects/def-ilafores/ronj2303` |
| `/project/def-ilafores/ronj2303/Horse/results` (`$RESULTS`) | Final outputs. Project limit is 1 TB; scratch has 20 TB per user but can be purged |

Paths are in `scripts/fir/fir_env.sh`. Source it in every session, because the job scripts need `ILAFORES` and `MC` in the environment.

## Branches
- `main` (`merged_tables` until its PR is merged): general changes.
- `fir`: created by `setup_fir.sh`, in both `MethodsComparison` and `ILL_pipelines`, for Fir-only changes (Horse in `datasets.tsv`, one hardcoded container path in the ILL MetaPhlAn script). These are left uncommitted for the user.

## Order of work
1. `setup_fir.sh` (fill in `SRC`, and `HORSE_INDEX` only if the index is built on iv12), in tmux.
2. Paste the blocks of `scripts/fir/horse.sh` one at a time: sample sheet, horse index, kneaddata, sample tsv, then tools.
3. `check_all Horse` between blocks.
4. Merge, then rsync to `$RESULTS`.

## Current state
Nothing has run on Fir. None of `scripts/fir/` has run on a cluster yet, so treat the first task of each tool as a pilot: check its log, output and memory before trusting the rest.

## Monitoring
- Start in the background: `bash $MC/scripts/fir/watch.sh <YYYY-MM-DD> <first job id>`. It exits on a new problem or when no job is left.
- On `PROBLEM`: read the log in `$MC/logs/` (kneaddata: `data/Horse/preproc/logs/`), diagnose, tell the user, then rerun `watch.sh`. Reported tasks are remembered in `~/watch_seen.txt`.
- `check_all Horse` cells are `ok/empty/unfinished`: green, yellow (some empty), red (some unfinished). `check_all -m Horse` lists unfinished rows for `--array`.
- `sacct -j <id> --format=JobID,State,ExitCode,Elapsed,TotalCPU,MaxRSS -P` for memory use.

## Rules from the user
- Don't commit or push unless asked. Never `git reset --hard` or `git clean`.
- Don't rerun profilers unless asked. Never delete per-sample profiler outputs.
- Light commands only on login nodes; heavy work goes through `sbatch`.
- Don't spend time on the Olive dataset, GTDB 214 or GenBank.
- Report outcomes plainly, including anything skipped or untested.

## Lessons from the HPC
- **`gather_SLURM_fast.sh`:** an empty `SIG` once made `sketch` overwrite a sample's `paired_1.fastq.gz`, on the first task of a new dataset. Fixed, but pilot it on the first task and check the reads are intact.
- **Sourmash timeouts:** old copies left cut-off `_gather.csv` files that looked finished. The script now writes `*_gather.tmp.csv` and renames on success. This fix was still waiting for its first test when this file was written.
- **mOTUs:** only the last `-s` counts, so unmatched reads are merged into one file. Output goes to `MOTUS` (v3) and `MOTUS4`.
- **Host removal:** kneaddata gets human and horse as two indexes. Several indexes need the `kneaddata_gzip` commit that added them (`--db "index1 index2"`); `setup_fir.sh` stops if the clone lacks it. Each index sees all reads, which is slower than `--serial` but simpler, and `<sample>_contaminants/` gets one set of files per host.
- **Kraken:** the index isn't reloaded per sample. The silent wait is classification, then a single-threaded report step of about 2 minutes on complex samples. More threads don't help that step.
- **`/dev/shm`:** counts against job memory, is private to a job, and was half of RAM on the HPC nodes. `kraken_node.sh` falls back to local NVMe.
- **`sbatch --export=ALL,VAR=x`:** such jobs were cancelled by root on the HPC. Pass script arguments.
- **`#SBATCH -D/-o`:** the tool scripts hardcode HPC paths, overridden by `run()` in `horse.sh`. Any new sbatch call needs the same.

## To confirm on Fir
- Per-user job limit and maximum array size. Each tool array is one job per sample; `CHUNK` splits them.
- `df -h /dev/shm`. Local disk is 7.84 TB.
- Whether `StdEnv/2020`, `apptainer/1.1.5`, `seqkit`, `python/3.8.10`, `mugqic/bwa` load (`motus3_SLURM.sh` needs `/cvmfs/soft.mugqic`). If the venv fails on Zen 5, mOTUs 3 may need rebuilding.
- That `kneaddata`'s Trimmomatic adapter path in `/cvmfs/soft.mugqic` is visible on compute nodes.
- That ILL's MetaPhlAn script binds scratch paths in its container.
- `kraken_node.sh` memory with 16 workers: raise `WORKERS` only if the first run leaves room. Each worker can reach about 10 GB.
