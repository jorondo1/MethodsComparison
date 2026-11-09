#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH -A def-ilafores

# Sketch one sample alone, so that sketching gets its own sacct record. Same commands as sylph_SLURM.sh and
# gather_SLURM_fast.sh, which then find the sketches and skip that step.
# Usage: sbatch --array=1-N -n <threads> --mem=<mem> benchmark_sketch.sh <sylph|sourmash> <tools tsv> <output>
set -eo pipefail

TOOL=$1
OUT=$3
IFS=$'\t' read -r SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2 <<< "$(awk "NR==$SLURM_ARRAY_TASK_ID" "$2")"
echo "$TOOL sketch of $SAM_ID on $(hostname)"

case $TOOL in
	sylph)
		export PATH=$ILAFORES/programs/conda_envs/sylph/bin:$PATH
		mkdir -p $OUT/sylph_sketches_c200
		sylph sketch -1 $FQ_P1 -2 $FQ_P2 -S $SAM_ID -c 200 -t $SLURM_NTASKS -d $OUT/sylph_sketches_c200
		;;
	sourmash)
		module load StdEnv/2020 apptainer/1.1.5
		mkdir -p $OUT/signatures
		singularity exec --writable-tmpfs -e -B $(dirname $FQ_P1):$(dirname $FQ_P1) -B $ILAFORES:$ILAFORES \
			$ILAFORES/programs/ILL_pipelines/containers/sourmash.4.8.11.sif \
			sourmash sketch dna -p k=31,scaled=1000,abund --merge $SAM_ID \
			-o $OUT/signatures/$SAM_ID.sig $FQ_P1 $FQ_P2 $FQ_U1 $FQ_U2
		;;
	*) echo "Unknown tool: $TOOL"; exit 1;;
esac
