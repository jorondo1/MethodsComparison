#!/bin/bash

#SBATCH -D /project/def-ilafores/ronj2303/MethodsComparison
#SBATCH -o /project/def-ilafores/ronj2303/MethodsComparison/logs/mOTU-%A_%a.slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=46G
#SBATCH -N 1
#SBATCH -n 192
#SBATCH -A def-ilafores
#SBATCH -J motus4

echo "initializing variables..."

export HOME_DIR=/project/def-ilafores/ronj2303
export OUT_DIR="${1}"/MOTUS4
export SAM_LIST="${2}"
export SAM_NUM=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${SAM_LIST})
IFS=$'\t' read -r SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2 <<< "$SAM_NUM" # array it up
export SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2

echo "Copying $SAM_ID files :\
	$FQ_P1 \
	$FQ_P2 \
	$FQ_U1 \
	$FQ_U2"

if [[ ! -f ${OUT_DIR}/${SAM_ID}_profile.txt ]]; then

# Sorts fasta and removes the /1 suffix from Kneaddata. mOTUs4 is pickier than v3...
# After seqkit processing, unload the module to free memory
ml seqkit
for var in FQ_P1 FQ_P2 FQ_U1 FQ_U2; do
    OUTFILE="${SLURM_TMPDIR}/$(basename ${!var} .gz)"
    seqkit sort -n -j 8 "${!var}" | \
    seqkit replace -p '\/[12]$' -r '' > "${OUTFILE}"
    export "$var"="${OUTFILE}"
done

export FQ_U=${SLURM_TMPDIR}/FQ_U.fastq
cat $FQ_U1 $FQ_U2 > $FQ_U

# Unload seqkit to free any cached memory
module unload seqkit

#  Debugging:
echo "File sizes before mOTUs:"
ls -lh $FQ_P1 $FQ_P2 $FQ_U
j
# Now run mOTUs with clean memory
echo "copying mOTUs container..."
cp $HOME_DIR/ILL_pipelines/containers/mOTUs_v4.0.4.sif $SLURM_TMPDIR

mkdir -p $OUT_DIR
echo "output will be stored in $OUT_DIR"
echo "executing mOTUs..."

module load StdEnv/2020 apptainer/1.1.5 

singularity exec --writable-tmpfs -e \
-B $SLURM_TMPDIR:$SLURM_TMPDIR \
-B $HOME_DIR:$HOME_DIR \
$SLURM_TMPDIR/mOTUs_v4.0.4.sif \
motus profile -f $FQ_P1 -r $FQ_P2 -s ${FQ_U} \
-n $SAM_ID -t $SLURM_NTASKS \
-o ${SLURM_TMPDIR}/${SAM_ID}_profile.txt

cp ${SLURM_TMPDIR}/${SAM_ID}_profile.txt ${OUT_DIR}/${SAM_ID}_profile.txt
fi