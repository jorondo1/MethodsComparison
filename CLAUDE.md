# MethodsComparison: instructions for Claude

This repo is worked on from two machines, each with its own Claude Code session, and the two sessions don't share memory:
- **HPC:** `/net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison` (also mounted as `/net/nfs-bio/jbod2/...`). Profilers, merging, and the heavy R runs (`scripts/0_run_all.R`) happen here.
- **Mac:** RStudio. Gets data and results from the HPC through `scripts/utils/sync_from_hpc.sh` (rsync), and makes the figures.

## At the start of every session
1. Read `doc/HANDOFF.md`: what the other machine did last, and what's next.
2. Read `doc/TODO.md`: the open items.
3. Run `git log --oneline -10` and `git status` to see whether the other side committed something since the last HANDOFF entry.

## Before ending a chunk of work (or when the user switches machines)
- Update `doc/HANDOFF.md`:
  - rewrite "Current state" so it's true now;
  - add a dated entry at the top of the log, tagged `[HPC]` or `[Mac]`;
  - keep it short and give file paths.
- Tick or add items in `doc/TODO.md`.
- `doc/` is gitignored, so stage these files with `git add -f doc/HANDOFF.md doc/TODO.md`. The user commits and pushes; don't do it without asking.

## Standing rules
- Never `git reset --hard` or `git clean`: `data/` and `Out/` hold untracked results. To sync clones, use `git fetch` + `git pull` on the working branch, or a mixed reset.
- Don't rerun the profilers unless asked.
- On the HPC, run only light commands on the login node; heavy work goes through `salloc` or `sbatch`. R: `module load StdEnv/2023 r/4.4.0`.
- Never delete per-sample profiler outputs.
- `data/` and `Out/` are not tracked: they move by rsync only.
- Stage with `git add -u` (+ `-f` for `doc/`), not `git add .`: `rademu_tmp/` (HPC, 140 GB) is untracked.
