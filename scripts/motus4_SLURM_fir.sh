#!/bin/bash
#SBATCH -D /project/def-ilafores/ronj2303/MethodsComparison
#SBATCH -o /project/def-ilafores/ronj2303/MethodsComparison/logs/mOTU-%A_%a.slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=46G
#SBATCH -N 1
#SBATCH -c 48
#SBATCH -A def-ilafores
#SBATCH -J motus4

echo "initializing variables..."
export HOME_DIR=/project/def-ilafores/ronj2303
export OUT_DIR="${1}"/MOTUS4
export SAM_LIST="${2}"
export SAM_NUM=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${SAM_LIST})
IFS=$'\t' read -r SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2 <<< "$SAM_NUM"
export SAM_ID

echo "Original file paths:"
echo "  FQ_P1: $FQ_P1"
echo "  FQ_P2: $FQ_P2"
echo "  FQ_U1: $FQ_U1"
echo "  FQ_U2: $FQ_U2"

if [[ ! -f ${OUT_DIR}/${SAM_ID}_profile.txt ]]; then

# Copy input files to local SLURM_TMPDIR first (much faster I/O)
echo "Copying input files to local node storage at $(date)..."
cp "$FQ_P1" "$SLURM_TMPDIR/" &
cp "$FQ_P2" "$SLURM_TMPDIR/" &
cp "$FQ_U1" "$SLURM_TMPDIR/" &
cp "$FQ_U2" "$SLURM_TMPDIR/" &
wait
echo "Copy complete at $(date)"

# Update variables to point to local copies
FQ_P1="$SLURM_TMPDIR/$(basename $FQ_P1)"
FQ_P2="$SLURM_TMPDIR/$(basename $FQ_P2)"
FQ_U1="$SLURM_TMPDIR/$(basename $FQ_U1)"
FQ_U2="$SLURM_TMPDIR/$(basename $FQ_U2)"

echo "Local file paths:"
echo "  FQ_P1: $FQ_P1"
echo "  FQ_P2: $FQ_P2"
echo "  FQ_U1: $FQ_U1"
echo "  FQ_U2: $FQ_U2"

ml seqkit

# Use 8 cores per seqkit job (1 CCD), process files sequentially
SEQKIT_THREADS=8

echo "Starting seqkit processing at $(date)"
for var in FQ_P1 FQ_P2 FQ_U1 FQ_U2; do
    OUTFILE="${SLURM_TMPDIR}/$(basename ${!var} .gz)"
    echo "Processing ${var}: $(basename ${!var}) at $(date)"
    time seqkit sort -n -j $SEQKIT_THREADS "${!var}" | \
    seqkit replace -p '\/[12]$' -r '' > "${OUTFILE}"
    export "$var"="${OUTFILE}"
done
echo "Seqkit processing complete at $(date)"

export FQ_U=${SLURM_TMPDIR}/FQ_U.fastq
cat $FQ_U1 $FQ_U2 > $FQ_U

module unload seqkit

echo "File sizes before mOTUs:"
ls -lh $FQ_P1 $FQ_P2 $FQ_U

echo "copying mOTUs container..."
cp $HOME_DIR/ILL_pipelines/containers/mOTUs_v4.0.4.sif $SLURM_TMPDIR
mkdir -p $OUT_DIR
echo "output will be stored in $OUT_DIR"

echo "executing mOTUs at $(date)..."
module load StdEnv/2020 apptainer/1.1.5 

singularity exec --writable-tmpfs -e \
-B $SLURM_TMPDIR:$SLURM_TMPDIR \
-B $HOME_DIR:$HOME_DIR \
$SLURM_TMPDIR/mOTUs_v4.0.4.sif \
motus profile -f $FQ_P1 -r $FQ_P2 -s ${FQ_U} \
-n $SAM_ID -t $SLURM_CPUS_PER_TASK \
-o ${SLURM_TMPDIR}/${SAM_ID}_profile.txt

echo "mOTUs complete at $(date)"
cp ${SLURM_TMPDIR}/${SAM_ID}_profile.txt ${OUT_DIR}/${SAM_ID}_profile.txt
fi