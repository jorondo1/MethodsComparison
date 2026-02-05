#!/bin/bash
#SBATCH -D /project/def-ilafores/ronj2303/MethodsComparison
#SBATCH -o /project/def-ilafores/ronj2303/MethodsComparison/logs/mOTU-%A_%a.slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=46G
#SBATCH -N 1
#SBATCH -c 48  # More realistic for this workload
#SBATCH -A def-ilafores
#SBATCH -J motus4

echo "initializing variables..."
export HOME_DIR=/project/def-ilafores/ronj2303
export OUT_DIR="${1}"/MOTUS4
export SAM_LIST="${2}"
export SAM_NUM=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${SAM_LIST})
IFS=$'\t' read -r SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2 <<< "$SAM_NUM"
export SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2

echo "Copying $SAM_ID files :\
	$FQ_P1 \
	$FQ_P2 \
	$FQ_U1 \
	$FQ_U2"

if [[ ! -f ${OUT_DIR}/${SAM_ID}_profile.txt ]]; then

# Process files in parallel using GNU parallel or background jobs
ml seqkit

# Option A: Process all 4 files simultaneously (recommended)
process_file() {
    local input="$1"
    local output="$2"
    seqkit sort -n -j $((SLURM_CPUS_PER_TASK / 4)) "${input}" | \
    seqkit replace -p '\/[12]$' -r '' > "${output}"
}

export -f process_file
export SLURM_CPUS_PER_TASK

OUTFILE_P1="${SLURM_TMPDIR}/$(basename ${FQ_P1} .gz)"
OUTFILE_P2="${SLURM_TMPDIR}/$(basename ${FQ_P2} .gz)"
OUTFILE_U1="${SLURM_TMPDIR}/$(basename ${FQ_U1} .gz)"
OUTFILE_U2="${SLURM_TMPDIR}/$(basename ${FQ_U2} .gz)"

# Process in parallel
process_file "$FQ_P1" "$OUTFILE_P1" &
process_file "$FQ_P2" "$OUTFILE_P2" &
process_file "$FQ_U1" "$OUTFILE_U1" &
process_file "$FQ_U2" "$OUTFILE_U2" &
wait

export FQ_P1="$OUTFILE_P1"
export FQ_P2="$OUTFILE_P2"
export FQ_U=${SLURM_TMPDIR}/FQ_U.fastq
cat "$OUTFILE_U1" "$OUTFILE_U2" > $FQ_U

module unload seqkit

echo "File sizes before mOTUs:"
ls -lh $FQ_P1 $FQ_P2 $FQ_U

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
-n $SAM_ID -t $SLURM_CPUS_PER_TASK \
-o ${SLURM_TMPDIR}/${SAM_ID}_profile.txt

cp ${SLURM_TMPDIR}/${SAM_ID}_profile.txt ${OUT_DIR}/${SAM_ID}_profile.txt
fi