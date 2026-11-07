# Horse on Fir: handoff for the Claude session there

You are on Fir (or another Alliance cluster) and have no memory from the HPC sessions. Read this first, then `CLAUDE.md`, `doc/TODO.md` and `doc/HANDOFF.md`. Your first job is to monitor the processing of the Horse dataset for the user.

## Goal
Run Horse (PRJNA1102860, feces metagenome, shallow shotgun) through kneaddata and all profilers, as an extra dataset for the methods comparison. The user supplies an accession list and picks 500–700 of the 2,778 runs. Median run: 1.05 Gbp. The project also has 30 `metagenome` and 30 `synthetic metagenome` runs, likely controls.

## Layout (all in `scripts/fir/fir_env.sh`; source it in every session)
| Variable | Path | Holds |
|---|---|---|
| `REPO` | `/project/6047190/ronj2303/MethodsComparison` | git repo, branch `fir`. `data/Horse/` receives only the profiler outputs, at the end |
| `MC` | `~/scratch/mc_work` | working copy: jobs run with `-D $MC` and write `data/Horse/` here; `scripts` is a symlink to `$REPO/scripts` |
| `ILL_PIPELINES` | `~/scratch/ILL_pipelines` | branch `kneaddata_gzip`; containers in `containers/` |
| `DBS` | `~/scratch/ref_dbs` | all databases, including `host_genomes/` and `metaphlan4_db/` |
| `RAW` | `~/scratch/Horse/raw` | raw reads |
Ignore `~/scratch/MethodsComparison`: it's an old clone from February.

## Branches
`fir` (in `REPO` only) holds the Fir changes: these scripts, Horse in `datasets.tsv`, and container paths taken from `$ILL_PIPELINES`. ILL needs no Fir branch.

## Order of work
1. `bash $REPO/scripts/fir/setup_fir.sh` in tmux. It pulls from hpc-bio over a single ssh connection.
2. Paste the blocks of `scripts/fir/horse.sh` one at a time: sample sheet, horse index (built from `~/scratch/host_genomes/*.fna.gz`), kneaddata, sample tsv, then the tools.
3. `check_all Horse` between blocks.
4. Merge, then rsync to `$REPO/data/Horse`.

## Raw reads (2026-11-07)
- **Selection:** `data/Horse/raw/sample_list.csv` has 772 IDs, matched against ENA's `sample_alias` (NCBI "Sample name").
  - 744 match, one WGS run each: 547 GB of fastq.gz.
  - 28 `*_2014` IDs only have tiny runs in another project (PRJNA1528879, probably coprocultures), so they're skipped.
- **Location:** `RAW` is now `~/scratch/Horse/raw` (`fir_env.sh`), because /project doesn't have room for it. The filereport is `$RAW/filereport_Horse.tsv`.
- **Scripts** (`scripts/fir/`):
  - `download_sra.sh <filereport> [P=4]` runs on the login node, in tmux session `horse_dl`. NCBI S3 gives about 40 MB/s per stream.
  - `convert_sra.sh` is the conversion array: `sbatch -A def-ilafores -c 4 --mem 8G --time 1:00:00 --array=1-744 -o logs/convert-%A_%a.out $F/convert_sra.sh filereport_Horse.tsv`, run from `$RAW`. Pilot 2829682: 56 s, 3.6 GB MaxRSS.
- **Check:** `cut -f2 $RAW/download_status.tsv | sort | uniq -c` should show 744 ok.

## Current state
Raw read download for Horse started 2026-11-07 (see above). Nothing else has run on Fir. None of `scripts/fir/` has run on a cluster yet, so treat the first task of each tool as a pilot: check its log, output and memory before trusting the rest.

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
- **Never replace a script that running jobs may be reading,** even with an atomic rename. On NFS, compute nodes can lose the old file: bash stops silently where it was, with exit 0, and the job shows `COMPLETED` with no output. 22 MetaPhlAn profiles were lost this way on 2026-11-06. Wait for the jobs to end, or use a new file name. Scripts passed directly to `sbatch` are spooled and safe; scripts they call are not. Compare outputs against tasks: `COMPLETED` alone is not proof.
- **Old nodes can hang at `module load`:** logs stop growing while the task stays `RUNNING`, sometimes for a day. Mammoth nodes are ageing and this is expected; `bio` and `c-3t` are newer.
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
