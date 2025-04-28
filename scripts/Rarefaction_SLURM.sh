#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH -D /net/nfs-ip34/home/def-ilafores/analysis/MethodsComparison
#SBATCH -o /net/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/logs/Rarefaction-%A_%a.slurm.out
#SBATCH -e /net/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/logs/Rarefaction-%A_%a.slurm.out
#SBATCH --time=6:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --mem=500G
#SBATCH -A def-ilafores
#SBATCH -J rarefaction

module load StdEnv/2023 r/4.4.0

echo "Begin..."
Rscript $ANCHOR$PWD/scripts/9.1_Rarefaction_curves.R -t 48 -c 6 \
-i $ANCHOR$PWD/Out/ps_full.ls.RDS \
-o $ANCHOR$PWD/Out/Rarefaction.RDS \
--steps 50 --repeats 10
echo "Done!"