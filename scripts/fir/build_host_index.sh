#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH -A def-ilafores
#SBATCH -J host_index

# Bowtie2 index of one host genome.
# Usage: sbatch build_host_index.sh <genome FASTA, gzipped or not> <output prefix>
set -eo pipefail

module load bowtie2
mkdir -p "$(dirname "$2")"
bowtie2-build --threads $SLURM_CPUS_PER_TASK --seed 1 "$1" "$2"
