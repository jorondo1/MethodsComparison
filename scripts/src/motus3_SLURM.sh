#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH -D /net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison
#SBATCH -o /net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison/logs/mOTU-%A_%a.slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=31G
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A def-ilafores
#SBATCH -J motus3

set -eo pipefail

# Load motus

# load all modules first, activate the venv last so its bin/ stays first in PATH
module use /cvmfs/soft.mugqic/CentOS6/modulefiles
module load StdEnv/2020 python/3.8.10 mugqic/bwa
source /net/nfs-bio/jbod2/def-ilafores/programs/motus3_env/bin/activate

# make sure motus comes from the venv, not from a pip --user install in ~/.local/bin
if [[ "$(command -v motus)" != "$VIRTUAL_ENV"/* ]]; then
    echo "motus not found in venv $VIRTUAL_ENV on $(hostname) (found: $(command -v motus || echo none))"
    exit 1
fi

export OUT_DIR=${PWD}/"${1}"/MOTUS
export SAM_LIST="${2}"
export SAM_NUM=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${SAM_LIST})
IFS=$'\t' read -r SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2 <<< "$SAM_NUM" # array it up
export SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2

echo "processing $SAM_ID files:\
	$FQ_P1 \
	$FQ_P2 \
	$FQ_U1 \
	$FQ_U2"

mkdir -p $OUT_DIR

motus profile -db $ILAFORES/ref_dbs/motus3/db_mOTU \
	-f $FQ_P1 -r $FQ_P2 -s $FQ_U1 -s $FQ_U2 -n $SAM_ID \
	-t $SLURM_NTASKS -c -p -k mOTU -q | awk -F'\t' 'NR==1 || $4 != 0' > ${OUT_DIR}/${SAM_ID}_profile.txt