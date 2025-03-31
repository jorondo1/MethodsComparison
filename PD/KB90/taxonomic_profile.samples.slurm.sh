#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH -D /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/PD/KB90
#SBATCH -o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/logs/taxonomic_profile-%A_%a.slurm.out
#SBATCH --time=72:00:00
#SBATCH --mem=120G
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A def-ilafores
#SBATCH -J taxonomic_profile


echo "loading env"
module load StdEnv/2020 apptainer/1.1.5

export __sample_line=$(cat /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/PD/preproc/preprocessed_reads.sample.tsv | awk "NR==$SLURM_ARRAY_TASK_ID")
export __sample=$(echo -e "$__sample_line" | cut -f1)
export __fastq_file1=$(echo -e "$__sample_line" | cut -f2)
export __fastq_file2=$(echo -e "$__sample_line" | cut -f3)

sleep $[ ( $RANDOM % 90 ) + 1 ]s

bash /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/ILL_pipelines/scripts/taxonomic_profile.sample.sh \
-o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/PD/KB90 \
-tmp $SLURM_TMPDIR \
-t 24 -m 120G \
-s $__sample \
-fq1 $__fastq_file1 \
-fq2 $__fastq_file2 \
--kraken_db /nfs3_ib/nfs-ip34/home/def-ilafores/ref_dbs/kraken2_dbs/k2_standard_20241228 \
--bracken_readlen 150 \
--confidence 0.90


