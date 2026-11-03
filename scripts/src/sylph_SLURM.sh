#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH -D /net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison
#SBATCH -o /net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison/logs/sylph-%A_%a.slurm.out
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -A def-ilafores
#SBATCH -J sylph

# Usage: sbatch --array=1-N sylph_SLURM.sh <dataset> <sample tsv> <sylph db (.syldb/.syl2db)> <sylph-tax taxonomy> <db name>
#   db name: output goes to data/<dataset>/SYLPH_<db name>
#   taxonomy: a custom taxonomy tsv(.gz) (genome file name <tab> d__;p__;...;s__) or a sylph-tax identifier (e.g. GTDB_r232)
# Writes to data/<dataset>/SYLPH_<db name>/per_sample/:
#   <sample>_profile.tsv : sylph profile, Sequence_abundance column = estimated read counts per genome
#   <sample>_profile.log : sylph profile log (classification rate)
#   <sample>.sylphmpa    : sylph-tax profile, sequence_abundance column = estimated read counts per clade
# Once all samples are done, build the merged tables with sylph_merge.sh.
start_time=$(date +%s)
set -eo pipefail

export SYLPH=${SYLPH:-sylph} # sylph 1.0.0, from $SYLPH_ENV (below)
export SCALE=200 # must not be larger than the -c used to build the database

# Parse options and variables
export SAM_LIST="${2}"
export SYLPH_DB="${3}"
export TAX="${4}"
export DB_NAME="${5}"
export OUT_DIR=${PWD}/data/"${1}"/SYLPH_${DB_NAME}
export SAMPLE_DIR=${OUT_DIR}/per_sample
export SKETCH_DIR=${PWD}/data/"${1}"/sylph_sketches_c${SCALE}
echo "Exporting to $OUT_DIR"
echo "using $SYLPH_DB with taxonomy $TAX"

# sylph + sylph-tax: conda env under $ILAFORES, so every node sees it (unlike ~/miniconda3, missing on the old nodes).
# Putting its bin/ first in PATH is enough for these two command-line tools: no conda needed in the job.
export SYLPH_ENV=${SYLPH_ENV:-$ILAFORES/programs/conda_envs/sylph}
export PATH=$SYLPH_ENV/bin:$PATH
# fail early (missing command)
for cmd in $SYLPH sylph-tax; do
	command -v $cmd > /dev/null || { echo "$cmd not found on $(hostname) (SYLPH_ENV=$SYLPH_ENV)"; exit 1; }
done
echo "$($SYLPH --version), sylph-tax $(sylph-tax --version 2>&1 | tail -1)"

mkdir -p $SAMPLE_DIR $SKETCH_DIR

# Parse samples
export SAM_NUM=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${SAM_LIST})
IFS=$'\t' read -r SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2 <<< "$SAM_NUM" # array it up
export SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2

# Sketch sample (paired reads only: sylph can't combine pairs and singletons into one sample)
# -S names the sketch and the sample column after the sample instead of the fastq path
export SKETCH=${SKETCH_DIR}/${SAM_ID}.paired.sylsp
if [[ ! -f $SKETCH ]]; then
	echo "Sketch $SAM_ID"
	$SYLPH sketch -1 $FQ_P1 -2 $FQ_P2 -S $SAM_ID -c $SCALE -t $SLURM_NTASKS -d $SKETCH_DIR
else
	echo "Sample sketch found. Skipping..."
fi

# Profile
# --estimate-read-counts puts read counts in Sequence_abundance (forces -u); Taxonomic_abundance stays a percentage
export PROFILE=${SAMPLE_DIR}/${SAM_ID}_profile.tsv
if [[ ! -f $PROFILE ]]; then
	echo "Profile $SAM_ID against $SYLPH_DB..."
	# the log holds "% of reads detected in database", used by sylph_merge.sh for the classification rate
	$SYLPH profile --estimate-read-counts $SYLPH_DB $SKETCH -t $SLURM_NTASKS -o ${PROFILE}.tmp 2>&1 \
		| tee ${SAMPLE_DIR}/${SAM_ID}_profile.log
	mv ${PROFILE}.tmp $PROFILE # so an interrupted run isn't mistaken for a finished profile
else
	echo "Profile found. Skipping..."
fi

# Taxonomy: sums each genome's Sequence_abundance (here read counts) up the lineage
if [[ ! -f ${SAMPLE_DIR}/${SAM_ID}.sylphmpa ]]; then
	echo "Annotate taxonomy..."
	sylph-tax taxprof $PROFILE -t $TAX -o ${SAMPLE_DIR}/
else
	echo "Taxonomic profile found. Skipping..."
fi

echo "Done ! Elapsed time:"
end_time=$(date +%s)
echo $((end_time - start_time))
