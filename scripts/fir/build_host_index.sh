#!/bin/bash
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH -A def-ilafores
#SBATCH -J host_index

# Combined bowtie2 index, so that kneaddata removes both hosts with a single --db.
# Usage: sbatch build_host_index.sh <index 1> <index 2> <output prefix>. Index = bowtie2 prefix.
set -eo pipefail

module load StdEnv/2020 apptainer/1.1.5
kneaddata="singularity exec -e -B $ILAFORES:$ILAFORES $ILL_PIPELINES/containers/kneaddata.0.12.0.sif"

mkdir -p "$(dirname "$3")"
$kneaddata bowtie2-inspect "$1" > $SLURM_TMPDIR/host.fasta
$kneaddata bowtie2-inspect "$2" >> $SLURM_TMPDIR/host.fasta
singularity exec -e -B $ILAFORES:$ILAFORES -B $SLURM_TMPDIR:$SLURM_TMPDIR $ILL_PIPELINES/containers/kneaddata.0.12.0.sif \
	bowtie2-build --threads $SLURM_CPUS_PER_TASK --seed 1 $SLURM_TMPDIR/host.fasta "$3"
