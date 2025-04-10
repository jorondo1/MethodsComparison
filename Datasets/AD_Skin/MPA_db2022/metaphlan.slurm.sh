#!/bin/bash -l

#SBATCH --mail-type=END,FAIL
#SBATCH -D /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/AD_Skin/MPA_db2022
#SBATCH -o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/logs/metaphlan-%A_%a.slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=30G
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A def-ilafores
#SBATCH -J metaphlan


echo "loading env"
module load StdEnv/2020 apptainer/1.1.5

export __sample_line=$(cat /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/AD_Skin/preproc/preprocessed_reads.sample.tsv | awk "NR==$SLURM_ARRAY_TASK_ID")
export __sample=$(echo -e "$__sample_line" | cut -f1)
export __fastq_file1=$(echo -e "$__sample_line" | cut -f2)
export __fastq_file2=$(echo -e "$__sample_line" | cut -f3)
export __fastq_file1_single=$(echo -e "$__sample_line" | cut -f4)
export __fastq_file2_single=$(echo -e "$__sample_line" | cut -f5)

bash -l /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/ILL_pipelines/scripts/taxonomic_abundance.metaphlan.sh \
-o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/AD_Skin/MPA_db2022/ \
-tmp $SLURM_TMPDIR \
-t 24 \
-s $__sample \
-fq1 $__fastq_file1 \
-fq2 $__fastq_file2 \
-fq1_single $__fastq_file1_single \
-fq2_single $__fastq_file2_single \
-db /nfs3_ib/nfs-ip34/fast/def-ilafores/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212


