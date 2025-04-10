#!/bin/bash


#SBATCH --mail-type=END,FAIL
#SBATCH -D /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/Feces_RA/preproc
#SBATCH -o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/Feces_RA/preproc/logs/preprocess.kneaddata-%A_%a.slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=125G
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A def-ilafores
#SBATCH -J preprocess



echo "loading env"
module load StdEnv/2020 apptainer/1.1.5

export __sample_line=$(cat /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/Feces_RA/raw/samples_to_process.tsv | awk "NR==$SLURM_ARRAY_TASK_ID")
export __sample=$(echo -e "$__sample_line" | cut -f1)
export __fastq_file1=$(echo -e "$__sample_line" | cut -f2)
export __fastq_file2=$(echo -e "$__sample_line" | cut -f3)

bash /nfs3_ib/nfs-ip34/home/def-ilafores/programs/ILL_pipelines/scripts/preprocess.kneaddata.sh \
-o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/Feces_RA/preproc \
-tmp $SLURM_TMPDIR \
-t 24 -m 125G \
-s $__sample -fq1 $__fastq_file1 -fq2 $__fastq_file2 \
--trimmomatic_options "SLIDINGWINDOW:4:30 MINLEN:100" \
--trimmomatic_adapters "ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10" \
--bowtie2_options "--very-sensitive-local" \
--db /nfs3_ib/nfs-ip34/fast/def-ilafores/host_genomes/GRCh38_index/grch38_1kgmaj

