#!/bin/bash

# Slurm records of the benchmark jobs. Reads the job ids in <output>/jobs.tsv, writes <output>/sacct.tsv
# Usage: bash benchmark_gather.sh <output>. Every job id of a label is gathered; 10_benchmark_plot.R keeps the last COMPLETED job per label x row.
OUT=${1:?Usage: benchmark_gather.sh <output>}

sacct -j "$(cut -f2 "$OUT/jobs.tsv" | paste -sd,)" -P \
	--format=JobID,JobName,State,ExitCode,ElapsedRaw,TotalCPU,MaxRSS,AllocCPUS,ReqMem,NodeList,Start,End \
	> "$OUT/sacct.tsv"

echo "$(($(wc -l < "$OUT/sacct.tsv") - 1)) records in $OUT/sacct.tsv"
