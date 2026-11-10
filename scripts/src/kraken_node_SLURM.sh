#!/bin/bash
#SBATCH -N 1
#SBATCH -A def-ilafores
#SBATCH -J kraken

# Kraken2 + Bracken on one node, for one index and the given confidences. Rerun to resume: finished samples are skipped.
# Local (ip34) version of scripts/fir/kraken_node.sh. /dev/shm is private to the job and counts against --mem,
# so request the index size plus about 10 GB per worker (the single-threaded report step).
# The index is copied once per job; samples are split across WORKERS parallel kraken_local.sh loops.
# Usage: sbatch -p <partition> -c <cpus> --mem=<mem> --time=<t> kraken_node_SLURM.sh <index> <folder suffix, "" or _GTDB> <sample tsv> <output dir> "<confidences, e.g. 10 45 90>"
# A kreport left by a killed run is reused as if finished: move any kreport without bracken output aside before rerunning.
INDEX=$1
SUFFIX=$2
TSV=$3
OUT=$4
CONFS=${5:-10 45 90}
WORKERS=${WORKERS:-8}
THREADS=$(( SLURM_CPUS_PER_TASK / WORKERS ))
(( THREADS < 1 )) && THREADS=1

MC=${MC:-$ILAFORES/analysis/MethodsComparison}
module load StdEnv/2020 apptainer/1.1.5

echo "$(date '+%F %T') $(hostname): index $INDEX, confidences $CONFS, $WORKERS workers x $THREADS threads"
cp -r "$INDEX" /dev/shm/ || { echo "Index copy failed"; exit 1; }
DB=/dev/shm/$(basename "$INDEX")
echo "$(date '+%F %T') index copied"

TMP=${SLURM_TMPDIR:-/tmp/kraken_$SLURM_JOB_ID}
mkdir -p $TMP
split -n r/$WORKERS -d "$TSV" $TMP/chunk.

for conf in $CONFS; do
	for chunk in $TMP/chunk.*; do
		[[ -s $chunk ]] || continue
		bash $MC/scripts/src/kraken_local.sh --tsv $chunk --confidence 0.$conf \
			--output $OUT/KB${conf}${SUFFIX} --kraken_db $DB --threads $THREADS &
	done
	wait
	echo "$(date '+%F %T') confidence 0.$conf done"
done

rm -rf "$DB" $TMP/chunk.*
