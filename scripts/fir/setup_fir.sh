#!/bin/bash

# One-off setup on Fir: working dir, ILL branch, containers, databases, environments. About 720 GB: run in tmux on a login node.
# Safe to rerun: rsync resumes, and finished steps are left alone.
#   bash /project/6047190/ronj2303/MethodsComparison/scripts/fir/setup_fir.sh
# Needs ssh from Fir to SRC. One connection is shared by all the rsync calls, so you authenticate once.
set -eo pipefail

SRC=ronj2303@hpc-bio.ccs.usherbrooke.ca # sees /net/nfs-bio

source "$(dirname "$(realpath "$0")")/fir_env.sh"
IV12_IL=/net/nfs-bio/jbod2/def-ilafores
IV12_FAST=/net/nfs-bio/fast/def-ilafores
mkdir -m 700 -p ~/.ssh
export RSYNC_RSH="ssh -o ControlMaster=auto -o ControlPath=$HOME/.ssh/cm-%r@%h:%p -o ControlPersist=2h"
pull() { rsync -avS --partial --no-g --no-p "$@"; }

# Working dir: jobs run with -D $MC and write data/<dataset>/ there; scripts come from the repo
mkdir -p $MC/data $MC/logs $DBS/{kraken2_dbs,sourmash_db,sylph_db,motus3,metaphlan4_db,host_genomes} $ILAFORES/programs/conda_envs
[[ -e $MC/scripts ]] || ln -s $REPO/scripts $MC/scripts

# ILL: the multi-index kneaddata lives on kneaddata_gzip
git -C $ILL_PIPELINES checkout kneaddata_gzip
grep -q db_binds $ILL_PIPELINES/scripts/preprocess.kneaddata.sh ||
	{ echo "ILL_pipelines has no multi-index kneaddata: git pull in $ILL_PIPELINES"; exit 1; }

# Containers -----------------------------------------------------------------
for sif in kneaddata.0.12.0 metaphlan.4.2.4 mOTUs_v4.0.4.mp2 kraken.2.1.2 sourmash.4.8.11; do
	pull $SRC:$IV12_IL/programs/ILL_pipelines/containers/$sif.sif $ILL_PIPELINES/containers/
done

# Databases ------------------------------------------------------------------
pull $SRC:$IV12_FAST/metaphlan4_db/{2022,2023,2025} $DBS/metaphlan4_db/
pull $SRC:$IV12_FAST/host_genomes/GRCh38_index $DBS/host_genomes/
pull $SRC:$IV12_IL/ref_dbs/motus3/db_mOTU $DBS/motus3/
pull $SRC:$IV12_IL/ref_dbs/sourmash_db/{RefSeq_20250528.k31.sbt.zip,RefSeq_20250528.lineages.csv,gtdb-rs220-reps-k31.zip,gtdb-rs220.lineages.reps.csv} $DBS/sourmash_db/
pull $SRC:$IV12_IL/ref_dbs/sylph_db/{RefSeq_20250528_c200.syldb,RefSeq_20250528_taxonomy.tsv.gz,gtdb-r220-c200-dbv1.syldb,gtdb_r220_metadata.tsv.gz} $DBS/sylph_db/

# Environments: their scripts hardcode the iv12 path --------------------------
pull $SRC:$IV12_IL/programs/conda_envs/sylph $ILAFORES/programs/conda_envs/
pull $SRC:$IV12_IL/programs/motus3_env $ILAFORES/programs/
for env in $ILAFORES/programs/conda_envs/sylph $ILAFORES/programs/motus3_env; do
	grep -rIl "$IV12_IL/programs" $env/bin | xargs -r sed -i "s#$IV12_IL/programs#$ILAFORES/programs#g"
done

# Checks ---------------------------------------------------------------------
df -h /dev/shm
module load StdEnv/2020 apptainer/1.1.5 && singularity --version
$ILAFORES/programs/conda_envs/sylph/bin/sylph --version
$ILAFORES/programs/conda_envs/sylph/bin/sylph-tax --version
module load python/3.8.10 && $ILAFORES/programs/motus3_env/bin/motus --version
echo "Everything but Kraken is ready: the other profilers can be launched now."

# Kraken last: about 400 GB, the longest transfer -------------------------------
pull $SRC:$IV12_IL/ref_dbs/kraken2_dbs/{k2_standard_20241228,k2_gtdb_genome_reps_20241109} $DBS/kraken2_dbs/
echo "Setup finished, Kraken included."
