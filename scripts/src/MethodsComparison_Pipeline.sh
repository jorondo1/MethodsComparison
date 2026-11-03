export MC=$ILAFORES/analysis/MethodsComparison
export MAPLE_DATA=/net/nfs-bio/jbod2/def-ilafores/analysis/mapleRhizosphere/WMS
#export ILL_PIPELINES=$ILAFORES/analysis/MethodsComparison/ILL_pipelines
cd $MC
source scripts/src/myFunctions.sh

# Create <dataset>_TSV and NUM_<dataset> variables
# Choose the one to work with :
dataset_variables "P19_Saliva" "$PR19/Saliva/preproc/preprocessed_reads.sample.tsv" # √
dataset_variables "P19_Gut" "$PR19/Feces/preproc/preprocessed_reads.sample.tsv" # √
dataset_variables "Moss" "$MOSS/preproc/preprocessed_reads.sample.tsv" # √
dataset_variables "NAFLD" "$MC/data/NAFLD/preproc/preprocessed_reads.sample.tsv" #
dataset_variables "AD_Skin" "$MC/data/AD_Skin/preproc/preprocessed_reads.sample.tsv" # √
dataset_variables "PD" "$MC/data/PD/preproc/preprocessed_reads.sample.tsv" # 47-655
dataset_variables "Bee" "$MC/data/Bee/preproc/preprocessed_reads.sample.tsv" 
# dataset_variables "Olive" "$MC/data/Olive/preproc/preprocessed_reads.sample.tsv" #XXX
dataset_variables "RA_Gut" "$MC/data/RA_Gut/preproc/preprocessed_reads.sample.tsv" 
dataset_variables "Maple" "$MAPLE_DATA/preproc/preprocessed_reads.sample.tsv" 

################################################################################################################
# QC #################
################################################################################################################

mkdir -p  $DATASET_PATH/raw

# Download from ENA tsv file with grep subset ------
samples=($(grep "WGS"  $DATASET_PATH/raw/filereport_*_tsv.txt | cut -f1  | grep -v 'run'))
cd /fast2/def-ilafores/Olive/raw
grep -E "$(IFS="|"; echo "${samples[*]}")" " $DATASET_PATH/raw/" | bash

# Create TSV
:> $DATASET_PATH/raw/samples_to_process.tsv
for sample in "${samples[@]}"; do
	fq1=$(find $ANCHOR/fast2/def-ilafores/$DATASET/raw -type f -name "${sample}_1.fastq.gz")
	fq2=$(find $ANCHOR/fast2/def-ilafores/$DATASET/raw -type f -name "${sample}_2.fastq.gz")
	echo -e "${sample}\t${fq1}\t${fq2}" >>  $DATASET_PATH/raw/samples_to_process.tsv
done

#Process samples
mkdir -p /fast2/def-ilafores/$DATASET/preproc/logs
bash $ANCHOR/$ILAFORES/programs/ILL_pipelines/generateslurm_preprocess.kneaddata.sh \
	--sample_tsv $ANCHOR/$DATASET_PATH/raw/samples_to_process.tsv \
	--out $ANCHOR/fast2/def-ilafores/$DATASET/preproc \
	--trimmomatic_options "SLIDINGWINDOW:4:20 MINLEN:50" \
	--db $FAST/host_genomes/GRCh38_index/grch38_1kgmaj \
	--slurm_mem 30G --slurm_threads 24

# Find arrays of missing samples:
missing_samples=$(grep -n -v -f <(ls  $DATASET_PATH/preproc/*/*_1.fastq.gz | awk -F'/' '{print $3}')  $DATASET_PATH/raw/samples_to_process.tsv | cut -f1 -d: | tr '\n' ','); echo $missing_samples

rm -r  $DATASET_PATH/preproc/.throttle
sbatch --array="$missing_samples" /nfs3_ib/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison/PD/preproc/preprocess.kneaddata.slurm.sh
# gzip all, except if gz already there
find $(dirname $TSV) -type f -name '*fastq' -print0 | xargs -0 -P 8 -I {} bash -c 'if [ ! -f "{}.gz" ]; then pigz -k "{}"; fi'

# Replace all .fastq by .fastq.gz in tsv
for tsv in $(find $ILAFORES/analysis/ -maxdepth 5 -readable -name "preprocessed_reads.sample.tsv" 2>/dev/null); do
sed -i 's|.fastq.gz|.fastq|g' $tsv
sed -i 's|.fastq|.fastq.gz|g' $tsv
sed -i 's|/nfs3_ib/nfs-ip34/home|/net/nfs-ip34/jbod2|g' $tsv
sed -i 's|/net/nfs-ip34/home|/net/nfs-ip34/jbod2|g' $tsv
done


# Local not preprocessed dataset (Maple)  -----------------
mapfile -t samples < <(tail -n +2 ${MAPLE_DATA}/raw/OMX753-WGS-13-05-2025_IndexAssignment.csv | grep WGS - | cut -d, -f2 | sort -u)

:> $MAPLE_DATA/raw/samples_to_process.tsv
for sample in "${samples[@]}"; do
	fq1=$(find $MAPLE -type f -name "${sample}_R1.fastq.gz")
	fq2=$(find $MAPLE -type f -name "${sample}_R2.fastq.gz")
	echo -e "${sample}\t${fq1}\t${fq2}" >> $MAPLE_DATA/raw/samples_to_process.tsv
done

# download acer genome
mkdir -p $FAST/host_genomes/_genome_fasta/AcerSaccharum
for chromosome in CM061734.1 CM061735.1 CM061736.1 CM061737.1 CM061738.1 CM061739.1 CM061740.1 CM061741.1 CM061742.1 CM061743.1 CM061744.1 CM061745.1 CM061746.1; do
	efetch -db nucleotide -id ${chromosome} -format fasta > $FAST/host_genomes/_genome_fasta/AcerSaccharum/${chromosome}.fasta
done
# index 
ml bowtie2
cat $FAST/host_genomes/_genome_fasta/AcerSaccharum/*.fasta >  $FAST/host_genomes/_genome_fasta/AcerSaccharum/full_AcerSaccharum.fasta
bowtie2-build --threads 48 --seed 1 $FAST/host_genomes/_genome_fasta/AcerSaccharum/full_AcerSaccharum.fasta $FAST/host_genomes/AcerSaccharum_index/AcerSaccharum


mkdir -p $MAPLE_DATA/preproc
bash $ILAFORES/programs/ILL_pipelines/generateslurm_preprocess.kneaddata.sh \
	--sample_tsv $MAPLE_DATA/raw/samples_to_process.tsv \
	--out $MAPLE_DATA/preproc \
	--trimmomatic_options "SLIDINGWINDOW:4:20 MINLEN:50" \
	--db $FAST/host_genomes/AcerSaccharum_index/AcerSaccharum \
	--slurm_mem 60G --slurm_threads 24

sbatch --partition bio,c-3t,cpubase --array=1-88 /net/nfs-bio/jbod2/def-ilafores/analysis/mapleRhizosphere/WMS/preproc/preprocess.kneaddata.slurm.sh

missing_samples=$(grep -n -v -f <(ls  $MAPLE_DATA/preproc/*/*R1.fastq.gz | awk -F'/' '{print $11}')  $MAPLE_DATA/raw/samples_to_process.tsv | cut -f1 -d: | tr '\n' ','); echo $missing_samples

sbatch --partition bio,c-3t,cpubase --array=1 /net/nfs-bio/jbod2/def-ilafores/analysis/mapleRhizosphere/WMS/preproc/preprocess.kneaddata.slurm.sh

# gzip output
sed -i 's|.fastq|.fastq.gz|g' preproc/preprocessed_reads.sample.tsv
rm preproc/*/*.gz
# with safety if interrupted
find "$MAPLE_DATA/preproc" -type f -name '*fastq' -print0 \
| xargs -0 -P 12 -I {} bash -c '
    if [ -f "{}.gz" ] && gzip -t "{}.gz" 2>/dev/null; then
        :                      # valid, skip
    else
        rm -f "{}.gz"          # remove partial/corrupt
        pigz "{}"
    fi
'

dataset_variables "Maple" "$MAPLE_DATA/preproc/preprocessed_reads.sample.tsv" 


################################################################################################################
# mOTUs4  ########
################################################################################################################

sbatch --array=1-"$N_SAMPLES" --partition bio,c-3t ${MC}/scripts/src/motus4_SLURM.sh ${DATASET_PATH} ${TSV}
sbatch --array=1-"$N_SAMPLES" --partition bio,c-3t ${MC}/scripts/src/motus3_SLURM.sh ${DATASET_PATH} ${TSV}

check_output 'MOTUS4'  $DATASET_PATH _profile.txt
 
# Rerun missing MOTUS
sbatch --array="$FOUND" ${MC}/scripts/motus4_SLURM_fir.sh ${DATASET_PATH} ${TSV}.fir
sbatch --array="$FOUND" --mem 300G -c 64 ${MC}/scripts/motus4_SLURM_fir.sh ${DATASET_PATH} ${TSV}.fir
#PD : jobid 22101634
# number of species in db
motus_db='/jbod2/def-ilafores/programs/motu-profiler_env/lib/python3.8/site-packages/motus/db_mOTU'
cat $motus_db/db_mOTU_taxonomy_ref-mOTUs.tsv | cut -f2 | sort -u | wc
cat $motus_db/db_mOTU_taxonomy_meta-mOTUs.tsv | cut -f1 | sort -u | wc

# Classification rate
mkdir -p Out/classification_rates
# -F'\t': taxonomy names contain spaces, so splitting on whitespace puts $4 in the wrong column
cd $MC && find ./data/*/MOTUS* -name '*_profile.txt' -exec awk -F'\t' '
FNR == 1 { assigned = 0; unassigned = 0 }
/^#/ { next }
/^unassigned/ { unassigned = $4; next }
{ assigned += $4 }
ENDFILE { 
    total = unassigned + assigned;
    if (total > 0) {
        rate = assigned / total; # classified fraction, like the other tools
        printf "%s\t%.2f\t%.2f\t%.5f\n", FILENAME, assigned, unassigned, rate        
    } else {
        printf "%s\tNA\tNA\tNA\n", FILENAME
    }
}
' {} + > Out/classification_rates/motus_classification_rate.tsv

################################################################################################################
# MetaPhlAn4 ###
################################################################################################################

cd $ILL_PIPELINES; git checkout kneaddata_gzip; cd $MC
metaphlan="bash $ILL_PIPELINES/generateslurm_taxonomic_abundance.metaphlan.4.2.4.sh \
	--slurm_log $MC/logs --slurm_walltime 48:00:00"

# Generate SLURM scripts https://github.com/jflucier/ILL_pipelines/blob/main/generateslurm_taxonomic_abundance.metaphlan.sh
$metaphlan --slurm_threads 8 --slurm_mem 60G --sample_tsv $TSV --db $FAST/metaphlan4_db/2022/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $DATASET_PATH/MPA_db2022
sbatch --array=1-"$N_SAMPLES" --partition bio,c-3t ${DATASET_PATH}/MPA_db2022/metaphlan.slurm.sh

$metaphlan --slurm_threads 8 --slurm_mem 60G --sample_tsv $TSV --db $FAST/metaphlan4_db/2023/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $DATASET_PATH/MPA_db2023
sbatch --array=1-"$N_SAMPLES" --partition bio,c-3t ${DATASET_PATH}/MPA_db2023/metaphlan.slurm.sh

$metaphlan --slurm_threads 8 --slurm_mem 80G --sample_tsv ${TSV} --db $FAST/metaphlan4_db/2025/mpa_vJan25_CHOCOPhlAnSGB_202503 --out ${DATASET_PATH}/MPA_db2025
sbatch --array=1-"$N_SAMPLES" --partition c-3t,bio ${DATASET_PATH}/MPA_db2025/metaphlan.slurm.sh

# Rerun missing MPA
database="MPA_db2025"
missing_MPA=$(grep -n -v -f <(ls $DATASET_PATH/$database/*/*_profile.txt | xargs -I {} basename {} | sed 's/_profile\.txt//') $MAPLE_DATA/preproc/preprocessed_reads.sample.tsv | cut -f1 -d: | tr '\n' ','); echo $missing_MPA
sbatch --array="$missing_MPA" ${DATASET_PATH}/MPA_db2025/metaphlan.slurm.sh
# Bee missing: jobid 71631

# number of species by database
bzcat $FAST/metaphlan4_db/2022/mpa_vOct22_CHOCOPhlAnSGB_202212_species.txt.bz2 | cut -f2 | awk -F'|' '{print $7}' | sort -u | grep "s__" | wc
bzcat $FAST/metaphlan4_db/2023/mpa_vJun23_CHOCOPhlAnSGB_202307_species.txt.bz2 | cut -f2 | awk -F'|' '{print $7}' | sort -u | grep "s__" | wc
bzcat $FAST/metaphlan4_db/2025/mpa_vJan25_CHOCOPhlAnSGB_202503_species.txt.bz2 | cut -f2 | awk -F'|' '{print $7}' | sort -u | grep "s__" | wc

# Check completion status
check_output 'MPA_db2025'  $DATASET_PATH _profile.txt

# Classified % = sum of the species rows (s__, not t__). With unclassified estimation (--unclassified_estimation before 4.2,
# default since; check the command in line 2 of the profile) abundances are scaled by the mapped fraction, so each rank
# sums to the classified %. The profiles no longer have their UNCLASSIFIED line; 100 - sum(s__) matched it (31.82, AD_Skin SRR20002363).
cd $MC && find ./data/*/MPA* -name '*_profile.txt' -exec awk -F'\t' '
FNR == 1 { classified = 0 }
!/^#/ && $1 ~ /\|s__/ && $1 !~ /\|t__/ { classified += $3 }
ENDFILE { printf "%s\t%.2f\t%.2f\t%.5f\n", FILENAME, classified, 100 - classified, classified / 100 }
' {} + > Out/classification_rates/mpa_classification_rate.tsv

# Remove bowtie indexes
rm */MPA_db*/*/*.bowtie2.txt
rm */*/.throttle -r

# Post-process to keep only species lines and headers
find "$MC"/data/${DATASET}/MPA_* -type f -name "*_profile.txt" | \
  parallel -j 24 'grep -E "s__|#" {} > {}.temp && mv {}.temp {}'
  
################################################################################################################
# Kraken/bracken ###
# on /fast2/ local #
################################################################################################################

# Copy fastqs to /fast2
nice -n10 ionice -c2 -n7 rclone copy $MAPLE_DATA/preproc /net/nfs-ip34/fast2/def-ilafores/ --transfers 16 --checkers 4 --modify-window 5s --fast-list --no-update-modtime --retries 3 --retries-sleep 30s --low-level-retries 1 -v -L --size-only --exclude "*contam*";

for dir in /jbod2/def-ilafores/analysis/projet_PROVID19/Feces/preproc /jbod2/def-ilafores/analysis/projet_PROVID19/Saliva/preproc/jbod2/def-ilafores/analysis/boreal_moss/preproc; do
nice -n10 ionice -c2 -n7 rclone copy $dir /fast2/def-ilafores/preproc --transfers 16 --checkers 4 --modify-window 5s --fast-list --no-update-modtime --retries 3 --retries-sleep 30s --low-level-retries 1 -v -L --size-only --exclude "*contam*";
done

# rearrange tsvs to point to new path and remove anchors
export TSV=/net/nfs-ip34/fast2/def-ilafores/preproc/preprocessed_reads.sample.tsv
sed -i "s|${DIR}|/net/nfs-ip34/fast2/def-ilafores/preproc|g" ${TSV}

cp -r $ILAFORES/ref_dbs/kraken2_dbs/k2_standard_20241228 /dev/shm/
cp -r $ILAFORES/ref_dbs/kraken2_dbs/k2_gtdb_genome_reps_20241109 /dev/shm/

k2_std=/dev/shm/k2_standard_20241228 
k2_gtdb=/dev/shm/k2_gtdb_genome_reps_20241109
k2_local=$MC/scripts/src/kraken_local.sh

# RUN KRAKEN ON /fast2/
ml apptainer

bash $k2_local --tsv ${TSV} --confidence 0.10 --output  $DATASET_PATH/KB10 --kraken_db $k2_std --threads 36
bash $k2_local --tsv ${TSV} --confidence 0.45 --output  $DATASET_PATH/KB45 --kraken_db $k2_std --threads 36
bash $k2_local --tsv ${TSV} --confidence 0.90 --output  $DATASET_PATH/KB90 --kraken_db $k2_std --threads 36
bash $k2_local --tsv ${TSV} --confidence 0.10 --output  $DATASET_PATH/KB10_GTDB --kraken_db $k2_gtdb --threads 36
bash $k2_local --tsv ${TSV} --confidence 0.45 --output  $DATASET_PATH/KB45_GTDB --kraken_db $k2_gtdb --threads 36
bash $k2_local --tsv ${TSV} --confidence 0.90 --output  $DATASET_PATH/KB90_GTDB --kraken_db $k2_gtdb --threads 36

# Check completion status
check_output 'KB10 KB45 KB90 KB10_GTDB KB45_GTDB KB90_GTDB' "$DATASET_PATH" _bracken_S.MPA.TXT

# Record % sequence classification
mkdir -p Out/classification_rates
cd "$MC" && find . -name '*.kreport' ! -name '*bracken*' -exec awk '
FNR==1 {unclassified=$2; next} 
FNR==2 {classified=$2; rate=classified/(unclassified+classified); 
         printf "%s\t%d\t%d\t%.5f\n", FILENAME, classified, unclassified, rate}
' {} + > Out/classification_rates/kraken_classification_rate.tsv

database="KB90"
missing_KB=$(grep -n -v -f <(ls  $DATASET_PATH/$database/*/*/*_bracken_S.MPA.TXT | awk -F'/' '{print $3}' | sed 's/_profile\.txt//')  $DATASET_PATH/preproc/preprocessed_reads.sample.tsv | cut -f1 -d: | tr '\n' ','); echo $missing_KB
sbatch --array="$missing_KB"  $DATASET_PATH/$database/taxonomic_profile.samples.slurm.sh  $DATASET_PATH "\$${dataset}_TSV"

# Once completely done, remove heavy files from kraken out 
cd data/$DATASET
rm -r \
KB*/*/*_taxonomy_nt \
KB*/*/*/*.bracken \
KB*/*/*bugs_list.MPA.TXT \
KB*/*/*/*_temp.MPA.TXT \
KB*/*/*.kreport \
KB*/*/*_kronagrams \
*/.throttle

cd $MC

################################################################################################################
# Sourmash gather ###
################################################################################################################

sbatch --mem=120G -n 24 --array=1-"$N_SAMPLES" --partition bio,c-3t,cpubase $MC/scripts/src/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "genbank-2022.03"
sbatch --mem=31G -n 24 --array=1-"$N_SAMPLES" --partition bio,c-3t,cpubase $MC/scripts/src/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "/net/nfs-ip34/jbod2/def-ilafores/ref_dbs/sourmash_db/gtdb-rs214-reps.k31.zip"
sbatch --mem=80G -n 16 --array=1-"$N_SAMPLES" --partition bio,c-3t,cpubase $MC/scripts/src/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "/net/nfs-ip34/jbod2/def-ilafores/ref_dbs/sourmash_db/gtdb-rs214-full.k31.zip"
sbatch --mem=31G -n 24 --array=1-"$N_SAMPLES" --partition bio,c-3t,cpubase $MC/scripts/src/gather_SLURM_fast.sh "$DATASET" "$TSV" "/net/nfs-ip34/jbod2/def-ilafores/ref_dbs/sourmash_db/gtdb-rs220-reps-k31.zip"
# sbatch --mem=31G -n 24 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "/jbod2/def-ilafores/analysis/boreal_moss/genome_sketches/gtdb-rs214-rep-MAGs.sbt.zip"
sbatch --mem=31G -n 1 --array=1-"$N_SAMPLES" --partition bio,c-3t,cpubase $MC/scripts/src/gather_SLURM_fast.sh "$DATASET" "$TSV" "$ILAFORES/ref_dbs/sourmash_db/RefSeq_20250528.k31.sbt.zip" # Here using 1 see gather_SLURM_fast.sh script; pending sourmash 4.9 update to remove 

# Check completion status
check_output 'RefSeq_20250528' $DATASET_PATH _gather.csv

# Extract Sourmash lineage subset
#for i in P19_Saliva P19_Gut PD AD_Skin Moss RA_Gut Bee Olive NAFLD; do
# i=RA_Gut
for j in $(find data/$DATASET -maxdepth 1 -type d -name 'SM_*'); do
	db_name=$(basename $j)
	db=$(echo "$db_name" | cut -d'_' -f2 | cut -d'-' -f1,2)
	ver=$(echo "$db_name" | cut -d'_' -f3)

	if [[ ! -f $j/${db_name}_lineages.csv ]]; then
        cat $j/*_gather.csv | cut -d, -f10 | tail -n+2 | awk '{print $1}' | sed 's/"//g; s/^[^_]*_\([^.]*\).*$/\1/' | grep -v name | sort -u | grep -Fhf - $ILAFORES/ref_dbs/sourmash_db/${db}*${ver}*.lineages.csv > $j/${db_name}_lineages.csv
	fi
done
#done

# Compute classification rates from sourmash gather output : 
cd $MC && find ./data/*/SM* -name '*_gather.csv' -exec awk -F',' '
  FNR == 1 {
      for (i=1; i<=NF; i++) {
      if ($i == "f_unique_weighted") {
        col = i;
        break;
      }
    }
    next;
  }
  { sum += $col }
  END { if (col) printf "%s\t%.5f\n", FILENAME, sum }
' {} \; > Out/classification_rates/sourmash_classification_rate.tsv

################################################################################################################
# Sylph #######
################################################################################################################

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
$sylph sketch -l /net/nfs-ip34/fast2/def-ilafores/refseq_genomes/20250528_genome_list_anchor.txt  -c 200 -t 32 -o $ILAFORES/ref_dbs/sylph_db/RefSeq_20250528_c200

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
	}' $ILAFORES/ref_dbs/sourmash_db/RefSeq_20250528.lineages.csv \
	/net/nfs-ip34/fast2/def-ilafores/refseq_genomes/20250528_genome_list_anchor.txt \
	2> $ILAFORES/ref_dbs/sylph_db/RefSeq_20250528_no_taxonomy.txt \
	| gzip > $ILAFORES/ref_dbs/sylph_db/RefSeq_20250528_taxonomy.tsv.gz # keep .gz in the name: sylph-tax decompresses based on it
# Sanity check: lines in genome list = lines in taxonomy + lines in no_taxonomy (no_taxonomy = the 14 refseq plasmid.* bundles)

# GTDB Database -----------------------
# Pre-built db (gtdb-r220-c200-dbv1.syldb) from the sylph docs; sylph-tax ships its taxonomy files.
# Downloads every sylph-tax taxonomy (incl. gtdb_r220_metadata.tsv.gz) and records the folder in sylph-tax's config.
conda activate $ILAFORES/programs/conda_envs/sylph
sylph-tax download --download-to $ILAFORES/ref_dbs/sylph_db/

# Profile ----------------
# Don't run on mp2
# args: dataset, sample tsv, sylph db, sylph-tax taxonomy, output name (data/<dataset>/SYLPH_<name>/per_sample)
sbatch --array=1-"$N_SAMPLES" \
--partition bio,c-3t --mem=50G -n 2 \
$MC/scripts/src/sylph_SLURM.sh \
"$DATASET" "$TSV" $ILAFORES/ref_dbs/sylph_db/RefSeq_20250528_c200.syldb $ILAFORES/ref_dbs/sylph_db/RefSeq_20250528_taxonomy.tsv.gz "RefSeq"

sbatch --array=1-"$N_SAMPLES" \
--partition bio,c-3t --mem=50G -n 2 \
$MC/scripts/src/sylph_SLURM.sh \
"$DATASET" "$TSV" $ILAFORES/ref_dbs/sylph_db/gtdb-r220-c200-dbv1.syldb $ILAFORES/ref_dbs/sylph_db/gtdb_r220_metadata.tsv.gz "gtdb_rs220"

# Test globdb because ~30 samples are completely empty with refseq
sbatch --array=1-"$N_SAMPLES" --partition bio,c-3t $MC/scripts/src/sylph_SLURM.sh \
"$DATASET" "$TSV" $ILAFORES/ref_dbs/sylph_db/globdb_r232_sylph_c200.syldb $ILAFORES/ref_dbs/sylph_db/globdb_r232_taxonomy_sylph.tsv.gz  "glob_r232"

# test fungi
sbatch --array=1-"$N_SAMPLES" --partition bio,c-3t $MC/scripts/src/sylph_SLURM.sh \
"$DATASET" "$TSV" $ILAFORES/ref_dbs/sylph_db/fungi-refseq-2025-10-11-c200-dbv1.syldb $ILAFORES/ref_dbs/sylph_db/fungi_refseq_2025-10-11_metadata.tsv.gz  "fungi_refseq"

check_output 'RefSeq' $DATASET_PATH _profile.tsv


# Merged tables ----------------
# Once all array tasks are done: counts.tsv, taxonomy.tsv, merged_species.tsv in data/<dataset>/SYLPH_<name>
# Lists unfinished samples (and doesn't merge) or samples with no genome detected. Run from $MC.
# Counts are corrected with each sample's classification rate (see sylph_merge.sh), read from per_sample/<sample>_profile.log

for db in RefSeq gtdb_rs220; do
	bash $MC/scripts/src/sylph_merge.sh "$DATASET" $db "$TSV"
done

# Classification rates (% of reads detected in database), all datasets and dbs, same layout as the other tools
cat $MC/data/*/SYLPH_*/classification_rate.tsv > Out/classification_rates/sylph_classification_rate.tsv





