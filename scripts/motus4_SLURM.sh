#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH -D /net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison
#SBATCH -o /net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison/logs/mOTU-%A_%a.slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=31G
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A def-ilafores
#SBATCH -J motus4

echo "initializing variables..."

export OUT_DIR=${PWD}/"${1}"/MOTUS4
export SAM_LIST=$ANCHOR/"${2}"
export SAM_NUM=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${SAM_LIST})
IFS=$'\t' read -r SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2 <<< "$SAM_NUM" # array it up
export SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2
for var in FQ_P1 FQ_P2 FQ_U1 FQ_U2; do
    if [[ "${!var}" != ${ANCHOR}/* && "${!var}" != /nfs3_ib/* ]]; then
        declare "$var"="${ANCHOR}/${!var}"
    fi
done

echo "$SAM_ID files will be processed:\
	$FQ_P1 \
	$FQ_P2 \
	$FQ_U1 \
	$FQ_U2"

echo "copying mOTUs container..."
cp $ANCHOR/jbod2/def-ilafores/programs/ILL_pipelines/containers/mOTUs_v4.0.4.sif $SLURM_TMPDIR

mkdir -p $OUT_DIR
echo "output will be stored in $OUT_DIR"
echo "executing mOTUs..."

module load StdEnv/2020 apptainer/1.1.5

singularity exec --writable-tmpfs -e \
-B $SLURM_TMPDIR:$SLURM_TMPDIR \
-B ${OUT_DIR}:${OUT_DIR} \
$SLURM_TMPDIR/mOTUs_v4.0.4.sif \
motus profile -f $FQ_P1 -r $FQ_P2 -s $FQ_U1 -s $FQ_U2 -n $SAM_ID \
-t $SLURM_NTASKS -o ${OUT_DIR}/${SAM_ID}_profile.txt
