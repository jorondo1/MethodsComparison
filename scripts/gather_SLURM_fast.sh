#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH -D /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison
#SBATCH -o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/logs/sourmash-%A_%a.slurm.out
#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A def-ilafores
#SBATCH -J sourmash
start_time=$(date +%s)

# Load motus
export sourmash="singularity exec --writable-tmpfs -e -B ${ANCHOR}/home:${ANCHOR}/home ${ILAFORES}/programs/ILL_pipelines/containers/sourmash.4.8.11.sif sourmash"
module load StdEnv/2020 apptainer/1.1.5

# Parse options
export SAM_LIST="${2}"
export SM_DB="${3}"
export OUT_DIR=${PWD}/"${1}"/SM_${SM_DB}

# Parse samples
export SAM_NUM=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${SAM_LIST})
IFS=$'\t' read -r SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2 <<< "$SAM_NUM" # array it up
export SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2
export SIG="${OUT_DIR}/../signatures/${SAM_ID}.sig"

mkdir -p $OUT_DIR/../signatures

if [[ ! -f $SIG ]]; then
	echo "Sketch metagenomes"
	$sourmash sketch dna -p k=31,scaled=1000,abund --merge ${SAM_ID} -o $SIG $FQ_P1 $FQ_P2 $FQ_U1 $FQ_U2
else
	echo "Metagenome sketches found. Skipping..."
fi

if [[ ! -f ${OUT_DIR}/${SAM_ID}_${DB_NAME}_gather.csv ]]; then
	echo "Gather against index database"
	$sourmash scripts fastgather $SIG $ILAFORES/ref_dbs/sourmash_db/${SM_DB}*k31.zip \
	-o ${OUT_DIR}/${SAM_ID}_${DB_NAME}_gather.csv \
	-c $SLURM_NTASKS
else
	echo "Gather output found. Skipping..."
fi

echo "Done ! Elapsed time:"
end_time=$(date +%s)
echo $((end_time - start_time))