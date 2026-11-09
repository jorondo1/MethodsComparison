#!/bin/bash
#SBATCH -A def-ilafores

# Kraken2 + Bracken benchmark for one index. /dev/shm is private to a job, so a job per sample would copy the index
# each time: this job copies it once, then runs each sample x confidence as its own srun step, one at a time.
# sacct gives one record per step, named <prefix><confidence>_<row>, for example KG45_7.
# Usage: sbatch -c <threads> --mem=<mem> benchmark_kraken.sh <index> <prefix> <folder suffix> <tools tsv> <rows> <output>
INDEX=$1
PREFIX=$2
SUFFIX=$3
TSV=$4
ROWS=$5
OUT=$6

module load StdEnv/2020 apptainer/1.1.5
SHM=/dev/shm/$(basename $INDEX)

srun -n 1 -c 1 --exact --mem=0 --job-name=${PREFIX}_copy cp -r $INDEX /dev/shm/ || { echo "Index copy failed"; exit 1; }

for conf in 10 45 90; do
	for row in $(seq 1 $ROWS); do
		sed -n "${row}p" $TSV > $SLURM_TMPDIR/row.tsv
		srun -n 1 -c $SLURM_CPUS_PER_TASK --exact --mem=0 --job-name=${PREFIX}${conf}_$row \
			bash $MC/scripts/src/kraken_local.sh --tsv $SLURM_TMPDIR/row.tsv --confidence 0.$conf \
			--output $OUT/KB${conf}${SUFFIX} --kraken_db $SHM --threads $SLURM_CPUS_PER_TASK \
			|| echo "Failed: ${PREFIX}${conf}_$row" # the other samples can still run
	done
done
