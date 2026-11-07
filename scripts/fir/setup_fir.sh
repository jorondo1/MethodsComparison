#!/bin/bash

# One-off setup on Fir: repos, containers, databases, environments. About 720 GB: run in tmux on a login node.
# Safe to rerun: rsync resumes, and finished steps are left alone.
# Clone main, or merged_tables until its PR is merged:
#   git clone --single-branch -b main https://github.com/jorondo1/MethodsComparison.git ~/scratch/MethodsComparison
#   bash ~/scratch/MethodsComparison/scripts/fir/setup_fir.sh
set -eo pipefail

SRC=""         # user@host that Fir can reach and that sees /net/nfs-bio
HORSE_INDEX="" # folder name in host_genomes on the HPC, if the horse index is built there. Blank: build it on Fir with horse.sh

[[ -n $SRC ]] || { echo "Fill in SRC at the top of this script first"; exit 1; }

source "$(dirname "$(realpath "$0")")/fir_env.sh"
IV12_IL=/net/nfs-bio/jbod2/def-ilafores
IV12_FAST=/net/nfs-bio/fast/def-ilafores
PULL="rsync -avS --partial --no-g --no-p"

mkdir -p $ILAFORES/programs/conda_envs $ILAFORES/ref_dbs/{kraken2_dbs,sourmash_db,sylph_db,motus3} \
	$FAST/metaphlan4_db $FAST/host_genomes $MC/logs

# Repos: Fir-specific changes live on a local fir branch ---------------------
# An older clone may already be here: never overwrite it
branch=$(git -C $MC branch --show-current)
if [[ ! $branch =~ ^(main|merged_tables|fir)$ ]] || { [[ $branch != fir ]] && [[ -n $(git -C $MC status --porcelain --untracked-files=no) ]]; }; then
	echo "$MC is on $branch or has local changes. Sort that out first, then check out main."
	exit 1
fi
[[ -d $ILL_PIPELINES ]] || git clone -b kneaddata_gzip https://github.com/jflucier/ILL_pipelines.git $ILL_PIPELINES
grep -q db_binds $ILL_PIPELINES/scripts/preprocess.kneaddata.sh ||
	{ echo "ILL_pipelines has no multi-index kneaddata: push kneaddata_gzip from iv12, then git pull in $ILL_PIPELINES"; exit 1; }
for repo in $MC $ILL_PIPELINES; do
	git -C $repo rev-parse --verify -q fir > /dev/null || git -C $repo checkout -b fir
done
git -C $MC checkout fir
git -C $ILL_PIPELINES checkout fir

# Containers -----------------------------------------------------------------
for sif in kneaddata.0.12.0 metaphlan.4.2.4 mOTUs_v4.0.4.mp2 kraken.2.1.2 sourmash.4.8.11; do
	$PULL $SRC:$IV12_IL/programs/ILL_pipelines/containers/$sif.sif $ILL_PIPELINES/containers/
done

# Databases ------------------------------------------------------------------
$PULL $SRC:$OLD_FAST/metaphlan4_db/{2022,2023,2025} $FAST/metaphlan4_db/
$PULL $SRC:$OLD_FAST/host_genomes/GRCh38_index $FAST/host_genomes/
[[ -z $HORSE_INDEX ]] || $PULL $SRC:$OLD_FAST/host_genomes/$HORSE_INDEX $FAST/host_genomes/
$PULL $SRC:$OLD_IL/ref_dbs/kraken2_dbs/{k2_standard_20241228,k2_gtdb_genome_reps_20241109} $ILAFORES/ref_dbs/kraken2_dbs/
$PULL $SRC:$OLD_IL/ref_dbs/motus3/db_mOTU $ILAFORES/ref_dbs/motus3/
$PULL $SRC:$OLD_IL/ref_dbs/sourmash_db/{RefSeq_20250528.k31.sbt.zip,RefSeq_20250528.lineages.csv,gtdb-rs220-reps-k31.zip,gtdb-rs220.lineages.reps.csv,gtdb-rs220-full.lineages.csv} $ILAFORES/ref_dbs/sourmash_db/
$PULL $SRC:$OLD_IL/ref_dbs/sylph_db/{RefSeq_20250528_c200.syldb,RefSeq_20250528_taxonomy.tsv.gz,gtdb-r220-c200-dbv1.syldb,gtdb_r220_metadata.tsv.gz} $ILAFORES/ref_dbs/sylph_db/

# Environments: their scripts hardcode the old path --------------------------
$PULL $SRC:$IV12_IL/programs/conda_envs/sylph $ILAFORES/programs/conda_envs/
$PULL $SRC:$IV12_IL/programs/motus3_env $ILAFORES/programs/
for env in $ILAFORES/programs/conda_envs/sylph $ILAFORES/programs/motus3_env; do
	grep -rIl "$IV12_IL/programs" $env/bin | xargs -r sed -i "s#$IV12_IL/programs#$ILAFORES/programs#g"
done

# Fir-only changes -----------------------------------------------------------
grep -q '^Horse' $MC/scripts/src/datasets.tsv ||
	printf 'Horse\t$ILAFORES/MethodsComparison/data/Horse/preproc/preprocessed_reads.sample.tsv\n' >> $MC/scripts/src/datasets.tsv

# Checks ---------------------------------------------------------------------
df -h /dev/shm
module load StdEnv/2020 apptainer/1.1.5 && singularity --version
$ILAFORES/programs/conda_envs/sylph/bin/sylph --version
$ILAFORES/programs/conda_envs/sylph/bin/sylph-tax --version
module load python/3.8.10 && $ILAFORES/programs/motus3_env/bin/motus --version
echo "Setup finished. Next: scripts/fir/horse.sh, and doc/HANDOFF_FIR.md for the Claude session on Fir."
