#!/bin/bash
#SBATCH -N 1
#SBATCH -c 192
#SBATCH --mem=0
#SBATCH --time=2-00:00:00
#SBATCH -A def-ilafores
#SBATCH -J kraken

# Kraken2 + Bracken on a whole node, for one index and 3 confidences. Rerun to resume: finished samples are skipped.
# The samples are shallow, so most of the time is the single-threaded report step: many small workers beat one large run.
# Each worker can reach about 10 GB in that step, on top of the index, so check memory before raising WORKERS.
# Usage: sbatch kraken_node.sh <index> <folder suffix, "" or _GTDB> <sample tsv> <output>
INDEX=$1
SUFFIX=$2
TSV=$3
OUT=$4
WORKERS=${WORKERS:-16}
THREADS=${THREADS:-8} # one CCD

module load StdEnv/2020 apptainer/1.1.5

# /dev/shm if the index fits, otherwise the node's NVMe. On NVMe the index is still shared, through the page cache.
index_gb=$(du -s --apparent-size -BG "$INDEX" | cut -f1 | tr -d G)
shm_gb=$(df --output=avail -BG /dev/shm | tail -1 | tr -d ' G')
DEST=$SLURM_TMPDIR
(( shm_gb > index_gb + 10 )) && DEST=/dev/shm
echo "Index ${index_gb}G, /dev/shm free ${shm_gb}G: using $DEST"
cp -r "$INDEX" "$DEST/" || { echo "Index copy failed"; exit 1; }
DB=$DEST/$(basename "$INDEX")

split -n r/$WORKERS -d "$TSV" $SLURM_TMPDIR/chunk.

for conf in 10 45 90; do
	for chunk in $SLURM_TMPDIR/chunk.*; do
		bash $MC/scripts/src/kraken_local.sh --tsv $chunk --confidence 0.$conf \
			--output $OUT/KB${conf}${SUFFIX} --kraken_db $DB --threads $THREADS &
	done
	wait
done
