#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH -A def-ilafores

# Bases per read file, to normalise CPU time by input size.
# Usage: sbatch benchmark_bases.sh <output>. Reads <output>/samples.tsv, writes <output>/bases.tsv
set -eo pipefail

OUT=$1
ml seqkit
cut -f3-6 "$OUT/samples.tsv" | tr '\t' '\n' | xargs seqkit stats -T -j $SLURM_NTASKS > "$OUT/bases.tsv"
