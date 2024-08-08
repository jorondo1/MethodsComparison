#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH -D /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison
#SBATCH -o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/logs/mOTU-%A_%a.slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=10G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -A def-ilafores
#SBATCH -J motus

# Load motus
source /nfs3_ib/nfs-ip34/home/def-ilafores/programs/motu-profiler_env/bin/activate
ml mugqic/bwa

export OUT_DIR=${PWD}/"${1}"/mOTU_abund
export SAM_LIST="${2}"
export SAM_NUM=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${SAM_LIST})
IFS=$'\t' read -r SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2 <<< "$SAM_NUM" # array it up
export SAM_ID FQ_P1 FQ_P2 FQ_U1 FQ_U2

mkdir -p $OUT_DIR

motus profile -f $FQ_P1 -r $FQ_P2 -s $FQ_U1 -s $FQ_U2 -n $SAM_ID \
	-t $SLURM_NTASKS -u -c -p  -q -k species | awk -F'\t' 'NR==1 || $4 != 0' > ${OUT_DIR}/${SAM_ID}_profile.txt
