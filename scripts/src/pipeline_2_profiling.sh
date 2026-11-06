# Profiling: launch each tool, check outputs, classification rates, merged tables.
# Copy-paste collection: running or sourcing this file does nothing. QC is in pipeline_1_QC.sh.
return 0 2>/dev/null; [[ $0 == "${BASH_SOURCE[0]}" ]] && { echo "Copy-paste collection: paste the commands, don't run this file"; exit 1; }

export MC=$ILAFORES/analysis/MethodsComparison
cd $MC && source scripts/src/myFunctions.sh && source scripts/src/classification_rates.sh

PARTS=bio,c-3t # + cpubase for light jobs (old nodes)
DB_MPA=$FAST/metaphlan4_db
DB_SM=$ILAFORES/ref_dbs/sourmash_db
DB_SYLPH=$ILAFORES/ref_dbs/sylph_db

use_dataset Maple # names: scripts/src/datasets.tsv

################################################################
# mOTUs  ########
################################################################

sbatch --array=1-"$N_SAMPLES" --partition $PARTS ${MC}/scripts/src/motus4_SLURM.sh ${DATASET_PATH} ${TSV}
sbatch --array=1-"$N_SAMPLES" --partition $PARTS ${MC}/scripts/src/motus3_SLURM.sh ${DATASET_PATH} ${TSV}

# Rerun missing mOTUs (big samples: --mem 300G -c 64)
check_output 'MOTUS3 MOTUS4' $DATASET_PATH _profile.txt
sbatch --array="$FOUND" --partition $PARTS ${MC}/scripts/src/motus4_SLURM.sh ${DATASET_PATH} ${TSV}

# number of species in db
motus_db='/jbod2/def-ilafores/programs/motu-profiler_env/lib/python3.8/site-packages/motus/db_mOTU'
cat $motus_db/db_mOTU_taxonomy_ref-mOTUs.tsv | cut -f2 | sort -u | wc
cat $motus_db/db_mOTU_taxonomy_meta-mOTUs.tsv | cut -f1 | sort -u | wc

# Classification rate
rate_motus $DATASET

################################################################
# MetaPhlAn4 ###
################################################################

cd $ILL_PIPELINES; git checkout kneaddata_gzip; cd $MC
metaphlan="bash $ILL_PIPELINES/generateslurm_taxonomic_abundance.metaphlan.4.2.4.sh \
	--slurm_log $MC/logs --slurm_walltime 48:00:00"

# Generate SLURM scripts https://github.com/jflucier/ILL_pipelines/blob/main/generateslurm_taxonomic_abundance.metaphlan.sh
$metaphlan --slurm_threads 24 --slurm_mem 31G --sample_tsv $TSV --db $DB_MPA/2022/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $DATASET_PATH/MPA_db2022
sbatch --array=1-"$N_SAMPLES" --partition cpubase ${DATASET_PATH}/MPA_db2022/metaphlan.slurm.sh

$metaphlan --slurm_threads 12 --slurm_mem 60G --sample_tsv $TSV --db $DB_MPA/2023/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $DATASET_PATH/MPA_db2023
sbatch --array=1-"$N_SAMPLES" --partition cpubase ${DATASET_PATH}/MPA_db2023/metaphlan.slurm.sh

$metaphlan --slurm_threads 12 --slurm_mem 60G --sample_tsv ${TSV} --db $DB_MPA/2025/mpa_vJan25_CHOCOPhlAnSGB_202503 --out ${DATASET_PATH}/MPA_db2025
sbatch --array=1-"$N_SAMPLES" --partition cpubase ${DATASET_PATH}/MPA_db2025/metaphlan.slurm.sh

# Rerun missing MPA
check_output 'MPA_db2022 MPA_db2023 MPA_db2025' $DATASET_PATH _profile.txt
sbatch --array="$FOUND" --partition $PARTS ${DATASET_PATH}/$database/metaphlan.slurm.sh

# number of species by database
bzcat $DB_MPA/2022/mpa_vOct22_CHOCOPhlAnSGB_202212_species.txt.bz2 | cut -f2 | awk -F'|' '{print $7}' | sort -u | grep "s__" | wc
bzcat $DB_MPA/2023/mpa_vJun23_CHOCOPhlAnSGB_202307_species.txt.bz2 | cut -f2 | awk -F'|' '{print $7}' | sort -u | grep "s__" | wc
bzcat $DB_MPA/2025/mpa_vJan25_CHOCOPhlAnSGB_202503_species.txt.bz2 | cut -f2 | awk -F'|' '{print $7}' | sort -u | grep "s__" | wc

# Classification rate
rate_mpa $DATASET

# ── DESTRUCTIVE: only after merge_tables.sh and rate_mpa ──
# Remove bowtie indexes
rm $DATASET_PATH/MPA_db*/*/*.bowtie2.txt
rm -r $DATASET_PATH/*/.throttle

# Post-process to keep only species lines and headers
find $DATASET_PATH/MPA_* -type f -name "*_profile.txt" | \
  parallel -j 24 'grep -E "s__|#" {} > {}.temp && mv {}.temp {}'
# ────

################################################################################################################
# Kraken/bracken ###
# on /fast2/ local #
################################################################################################################

# On the allocated node (/dev/shm there is big enough)
cp -r $ILAFORES/ref_dbs/kraken2_dbs/k2_standard_20241228 /dev/shm/
cp -r $ILAFORES/ref_dbs/kraken2_dbs/k2_gtdb_genome_reps_20241109 /dev/shm/

k2_std=/dev/shm/k2_standard_20241228
k2_gtdb=/dev/shm/k2_gtdb_genome_reps_20241109
k2_local=$MC/scripts/src/kraken_local.sh

# RUN KRAKEN (finished samples are skipped, so rerun the same loop for missing ones)
ml apptainer
for conf in 10 45 90; do
	#bash $k2_local --tsv ${TSV} --confidence 0.$conf --output $DATASET_PATH/KB$conf --kraken_db $k2_std --threads $SLURM_CPUS_PER_TASK
	bash $k2_local --tsv ${TSV} --confidence 0.$conf --output $DATASET_PATH/KB${conf}_GTDB --kraken_db $k2_gtdb --threads $SLURM_CPUS_PER_TASK
done

# Check completion status
check_output 'KB10 KB45 KB90 KB10_GTDB KB45_GTDB KB90_GTDB' "$DATASET_PATH" _bracken_S.MPA.TXT

# Classification rate: needs the .kreport files, so before the cleanup below
rate_kraken $DATASET

# ── DESTRUCTIVE: only after merge_tables.sh and rate_kraken ──
# Remove heavy files from kraken out
cd $DATASET_PATH
rm -r \
KB*/*/*_taxonomy_nt \
KB*/*/*/*.bracken \
KB*/*/*bugs_list.MPA.TXT \
KB*/*/*/*_temp.MPA.TXT \
KB*/*/*_kronagrams \
*/.throttle
cd $MC
# ────

################################################################################################################
# Sourmash gather ###
################################################################################################################

sbatch --mem=120G -n 24 --array=1-"$N_SAMPLES" --partition $PARTS,cpubase $MC/scripts/src/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "genbank-2022.03"
sbatch --mem=31G -n 24 --array=1-"$N_SAMPLES" --partition $PARTS,cpubase $MC/scripts/src/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "$DB_SM/gtdb-rs214-reps.k31.zip"
sbatch --mem=80G -n 16 --array=1-"$N_SAMPLES" --partition cpubase $MC/scripts/src/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "$DB_SM/gtdb-rs214-full.k31.zip"
sbatch --mem=31G -n 24 --array=1-"$N_SAMPLES" --partition cpubase $MC/scripts/src/gather_SLURM_fast.sh "$DATASET" "$TSV" "$DB_SM/gtdb-rs220-reps-k31.zip"
# sbatch --mem=31G -n 24 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "/jbod2/def-ilafores/analysis/boreal_moss/genome_sketches/gtdb-rs214-rep-MAGs.sbt.zip"
sbatch --mem=31G -n 1 --array=1-"$N_SAMPLES" --partition $PARTS,cpubase $MC/scripts/src/gather_SLURM_fast.sh "$DATASET" "$TSV" "$DB_SM/RefSeq_20250528.k31.sbt.zip" # Here using 1 see gather_SLURM_fast.sh script; pending sourmash 4.9 update to remove

# Check completion status
check_output 'RefSeq_20250528' $DATASET_PATH _gather.csv

# Extract Sourmash lineage subset (only for the old per-sample R parsing; merge_tables.sh uses the full DB lineages)
for j in $(find $DATASET_PATH -maxdepth 1 -type d -name 'SM_*'); do
	db_name=$(basename $j)
	db=$(echo "$db_name" | cut -d'_' -f2 | cut -d'-' -f1,2)
	ver=$(echo "$db_name" | cut -d'_' -f3)

	if [[ ! -f $j/${db_name}_lineages.csv ]]; then
        cat $j/*_gather.csv | cut -d, -f10 | tail -n+2 | awk '{print $1}' | sed 's/"//g; s/^[^_]*_\([^.]*\).*$/\1/' | grep -v name | sort -u | grep -Fhf - $DB_SM/${db}*${ver}*.lineages.csv > $j/${db_name}_lineages.csv
	fi
done

# Classification rate
rate_sourmash $DATASET

################################################################################################################
# Sylph #######
################################################################################################################

# ── ONE-TIME (done): sylph env and databases ──
# sylph + sylph-tax env (once). Under $ILAFORES, not ~/miniconda3/envs: the old nodes don't see /home/ronj2303,
# and a conda env can't be moved afterwards (its scripts hardcode its path). The SLURM scripts only put its bin/ in PATH.
conda create -p $ILAFORES/programs/conda_envs/sylph -c conda-forge -c bioconda sylph=1.0.0 sylph-tax
$ILAFORES/programs/conda_envs/sylph/bin/sylph --version

# REFSEQ Database -----------------------
# [patch] Append anchor to genome files
cat /net/nfs-ip34/fast2/def-ilafores/refseq_genomes/20250528_genome_list.txt | \
sed 's|/fast2|/net/nfs-ip34/fast2|g' > \
/net/nfs-ip34/fast2/def-ilafores/refseq_genomes/20250528_genome_list_anchor.txt

# Sketch refseq custom db
$sylph sketch -l /net/nfs-ip34/fast2/def-ilafores/refseq_genomes/20250528_genome_list_anchor.txt  -c 200 -t 32 -o $DB_SYLPH/RefSeq_20250528_c200

# sylph-tax custom taxonomy for the refseq db: genome file name <tab> d__;p__;c__;o__;f__;g__;s__
# Lineages come from the sourmash lineage csv (accession,taxid,7 ranks; see sourmash_refseq_db.sh).
# A genome is matched on its accession, else on its 9-digit number (lineage listed under the GCA/GCF twin or another version).
# Genomes without a lineage go to RefSeq_20250528_no_taxonomy.txt; sylph-tax reports them as NO_TAXONOMY.
awk -F',' '
	NR == FNR {
		gsub(/"/, "")
		lineage = "d__" $(NF-6) ";p__" $(NF-5) ";c__" $(NF-4) ";o__" $(NF-3) ";f__" $(NF-2) ";g__" $(NF-1) ";s__" $NF
		tax[$1] = lineage
		split($1, a, "_"); by_num[substr(a[2], 1, 9)] = lineage
		next
	}
	{
		n = split($0, p, "/"); f = p[n]
		split(f, a, "_"); acc = a[1] "_" a[2]; num = substr(a[2], 1, 9)
		if (acc in tax) print f "\t" tax[acc]
		else if (num in by_num) print f "\t" by_num[num]
		else print f > "/dev/stderr"
	}' $DB_SM/RefSeq_20250528.lineages.csv \
	/net/nfs-ip34/fast2/def-ilafores/refseq_genomes/20250528_genome_list_anchor.txt \
	2> $DB_SYLPH/RefSeq_20250528_no_taxonomy.txt \
	| gzip > $DB_SYLPH/RefSeq_20250528_taxonomy.tsv.gz # keep .gz in the name: sylph-tax decompresses based on it
# Sanity check: lines in genome list = lines in taxonomy + lines in no_taxonomy (no_taxonomy = the 14 refseq plasmid.* bundles)

# GTDB Database -----------------------
# Pre-built db (gtdb-r220-c200-dbv1.syldb) from the sylph docs; sylph-tax ships its taxonomy files.
# Downloads every sylph-tax taxonomy (incl. gtdb_r220_metadata.tsv.gz) and records the folder in sylph-tax's config.
conda activate $ILAFORES/programs/conda_envs/sylph
sylph-tax download --download-to $DB_SYLPH/
# ────

# Profile ----------------
# Don't run on mp2
# args: dataset, sample tsv, sylph db, sylph-tax taxonomy, output name (data/<dataset>/SYLPH_<name>/per_sample)
sbatch --array=1-"$N_SAMPLES" \
--partition $PARTS --mem=50G -n 2 \
$MC/scripts/src/sylph_SLURM.sh \
"$DATASET" "$TSV" $DB_SYLPH/RefSeq_20250528_c200.syldb $DB_SYLPH/RefSeq_20250528_taxonomy.tsv.gz "RefSeq"

sbatch --array=1-"$N_SAMPLES" \
--partition $PARTS --mem=50G -n 2 \
$MC/scripts/src/sylph_SLURM.sh \
"$DATASET" "$TSV" $DB_SYLPH/gtdb-r220-c200-dbv1.syldb $DB_SYLPH/gtdb_r220_metadata.tsv.gz "gtdb_rs220"

# Test globdb because ~30 samples are completely empty with refseq
sbatch --array=1-"$N_SAMPLES" --partition $PARTS $MC/scripts/src/sylph_SLURM.sh \
"$DATASET" "$TSV" $DB_SYLPH/globdb_r232_sylph_c200.syldb $DB_SYLPH/globdb_r232_taxonomy_sylph.tsv.gz  "glob_r232"

# test fungi
sbatch --array=1-"$N_SAMPLES" --partition $PARTS $MC/scripts/src/sylph_SLURM.sh \
"$DATASET" "$TSV" $DB_SYLPH/fungi-refseq-2025-10-11-c200-dbv1.syldb $DB_SYLPH/fungi_refseq_2025-10-11_metadata.tsv.gz  "fungi_refseq"

check_output 'RefSeq' $DATASET_PATH _profile.tsv


# Merged tables ----------------
# Once all array tasks are done: counts.tsv, taxonomy.tsv, merged_species.tsv in data/<dataset>/SYLPH_<name>
# Lists unfinished samples (and doesn't merge) or samples with no genome detected. Run from $MC.
# Counts are corrected with each sample's classification rate (see sylph_merge.sh), read from per_sample/<sample>_profile.log

for db in RefSeq gtdb_rs220; do
	bash $MC/scripts/src/sylph_merge.sh "$DATASET" $db "$TSV"
done
# (or merge_tables.sh --tools sylph "$DATASET", which runs sylph_merge.sh on every SYLPH_* dir, see below)

# Classification rates (% of reads detected in database), from the classification_rate.tsv written by sylph_merge.sh
rate_sylph $DATASET

################################################################################################################
# Merged tables (all tools) #######
################################################################################################################

# One table per dataset x database, merged with each tool's own utility, read by scripts/0_proc_data.R
# (scripts/utils/parse_merged.R). Datasets and sample tsvs: scripts/src/datasets.tsv. Run from $MC.
# Without -f, folders with samples of the tsv missing are skipped (listed). See the header of merge_tables.sh for outputs.
sbatch $MC/scripts/src/merge_tables.sh -f                          # every dataset, every tool
bash $MC/scripts/src/merge_tables.sh --tools sourmash,kraken PD    # some tools, one dataset
