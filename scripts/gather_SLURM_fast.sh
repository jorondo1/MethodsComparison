#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH -D /net/nfs-ip34/home/def-ilafores/analysis/MethodsComparison
#SBATCH -o /net/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/logs/sourmash-%A_%a.slurm.out
#SBATCH --time=6:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A def-ilafores
#SBATCH -J sourmash
start_time=$(date +%s)

export ILAFORES=$ANCHOR/$ILAFORES
# Load sourmash
module load StdEnv/2020 apptainer/1.1.5

# Parse options and variables
export SAM_LIST="${ANCHOR}/${2}"
export SM_DB="${ANCHOR}/${3}"
export DB_NAME=$(basename ${SM_DB/[-._]k??.*/})
export OUT_DIR=${PWD}/"${1}"/SM_${DB_NAME}
export DB_ANCHOR=$(dirname $SM_DB)
echo "Exporting to $OUT_DIR"

# Parse samples
export SAM_NUM=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${SAM_LIST})
IFS=$'\t' read -r SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2 <<< "$SAM_NUM" # array it up
export SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2
# Add Anchor to sample paths
for var in FQ_P1 FQ_P2 FQ_U1 FQ_U2; do
    if [[ "${!var}" != ${ANCHOR}/* && "${!var}" != /nfs3_ib/* ]]; then
        declare "$var"="${ANCHOR}${!var}"
    fi
done

SAM_ANCHOR=$(dirname $FQ_P1)

export sourmash="singularity exec --writable-tmpfs -e -B ${SAM_ANCHOR}:${SAM_ANCHOR} -B $ILAFORES:$ILAFORES -B $DB_ANCHOR:DB_ANCHOR ${ILAFORES}/programs/ILL_pipelines/containers/sourmash.4.8.11.sif sourmash"

export SIG=$(realpath "${PWD}/${1}/signatures/${SAM_ID}.sig")

mkdir -p $OUT_DIR/../signatures

if [[ ! -f $SIG ]]; then
	echo "$SIG not found"
	echo "Sketch metagenomes"
	$sourmash sketch dna -p k=31,scaled=1000,abund --merge ${SAM_ID} -o $SIG $FQ_P1 $FQ_P2 $FQ_U1 $FQ_U2
else
	echo "Metagenome sketches found. Skipping..."
fi

if [[ ! -f ${OUT_DIR}/${SAM_ID}_${DB_NAME}_gather.csv ]]; then
	echo "Gather against index database"
	$sourmash scripts fastgather $SIG ${SM_DB} \
	-o ${OUT_DIR}/${SAM_ID}_${DB_NAME}_gather.csv \
	-c $SLURM_NTASKS
else
	echo "Gather output found. Skipping..."
fi

echo "Done ! Elapsed time:"
end_time=$(date +%s)
echo $((end_time - start_time))
