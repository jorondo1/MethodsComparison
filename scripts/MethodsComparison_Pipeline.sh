export MC=$ILAFORES/analysis/MethodsComparison
export ILL_PIPELINES=$ILAFORES/analysis/MethodsComparison/ILL_pipelines
cd $MC
source scripts/myFunctions.sh

# Create <dataset>_TSV and NUM_<dataset> variables
# Choose the one to work with :
dataset_variables "P19_Saliva" "$PR19/Saliva/preproc/preprocessed_reads.sample.tsv"
dataset_variables "P19_Gut" "$PR19/Feces/preproc/preprocessed_reads.sample.tsv"
dataset_variables "Moss" "$MOSS/preproc/preprocessed_reads.sample.tsv"
dataset_variables "NAFLD" "$MC/NAFLD/preproc/preprocessed_reads.sample.tsv"
dataset_variables "AD_Skin" "$MC/NAFLD/preproc/preprocessed_reads.sample.tsv"
dataset_variables "PD" "$MC/PD/preproc/preprocessed_reads.sample.tsv"
dataset_variables "Bee" "$MC/Bee/preproc/preprocessed_reads.sample.tsv"
dataset_variables "Olive" "$MC/Olive/preproc/preprocessed_reads.sample.tsv"
dataset_variables "RA_Gut" "$MC/RA_Gut/preproc/preprocessed_reads.sample.tsv"

######################
# QC #################
######################

mkdir -p $MC/$DATASET/raw

# Download from ENA tsv file with grep subset
samples=($(grep "WGS" $MC/$DATASET/raw/filereport_*_tsv.txt | cut -f1  | grep -v 'run'))
cd /fast2/def-ilafores/Olive/raw
grep -E "$(IFS="|"; echo "${samples[*]}")" "$MC/$DATASET/raw/" | bash

# Create TSV
:> $MC/$DATASET/raw/samples_to_process.tsv
for sample in "${samples[@]}"; do
	fq1=$(find $ANCHOR/fast2/def-ilafores/$DATASET/raw -type f -name "${sample}_1.fastq.gz")
	fq2=$(find $ANCHOR/fast2/def-ilafores/$DATASET/raw -type f -name "${sample}_2.fastq.gz")
	echo -e "${sample}\t${fq1}\t${fq2}" >> $MC/$DATASET/raw/samples_to_process.tsv
done

#Process samples
mkdir -p /fast2/def-ilafores/$DATASET/preproc/logs
bash $ANCHOR/$ILAFORES/programs/ILL_pipelines/generateslurm_preprocess.kneaddata.sh \
	--sample_tsv $ANCHOR/$MC/$DATASET/raw/samples_to_process.tsv \
	--out $ANCHOR/fast2/def-ilafores/$DATASET/preproc \
	--trimmomatic_options "SLIDINGWINDOW:4:20 MINLEN:50" \
	--db $FAST/host_genomes/GRCh38_index/grch38_1kgmaj \
	--slurm_mem 30G --slurm_threads 24
# correct script as the shell command needs to be anchored!
cp -r /fast2/def-ilafores/$DATASET/preproc/* $DATASET/preproc/

# Find arrays of missing samples:
missing_samples=$(grep -n -v -f <(ls $DATASET/preproc/*/*_1.fastq.gz | awk -F'/' '{print $3}') $DATASET/raw/samples_to_process.tsv | cut -f1 -d: | tr '\n' ','); echo $missing_samples

rm -r $DATASET/preproc/.throttle
sbatch --array="$missing_samples" /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/PD/preproc/preprocess.kneaddata.slurm.sh

################
# mOTUs ########
################

# Custom SLURM script
sbatch --array=1-"$N_SAMPLES" $MC/scripts/motus_SLURM.sh $DATASET $TSV.fast

# Check completion status
check_output 'MOTUS' $DATASET _profile.txt
 
# Rerun missing MOTUS
rm $DATASET/MOTUS/_profile.txt # not sure why that appears
missing_motus=$(grep -n -v -f <(ls "$DATASET/MOTUS/"*_profile.txt | awk -F'/' '{print $3}' | sed 's/_profile\.txt//') "$(eval echo \$${dataset}_TSV)" | cut -f1 -d: | tr '\n' ','); echo "$missing_motus"
sbatch --array="$missing_motus" $MC/scripts/motus_SLURM.sh $DATASET "$(eval echo \$${dataset}_TSV)"

####################
# Kraken/bracken ###
# on /fast2/ local #
####################

# Copy fastqs to /fast2
for dir in $(find $MC -maxdepth 3 -type d -name 'preproc'); do 
nice -n10 ionice -c2 -n7 rclone copy $dir /fast2/def-ilafores/preproc --transfers 16 --checkers 4 --modify-window 5s --fast-list --no-update-modtime --retries 3 --retries-sleep 30s --low-level-retries 1 -v -L --size-only --exclude "*contam*";
done

# rearrange tsvs to point to new path
for tsv in $(find $ILAFORES/analysis/ -maxdepth 4 -name "preprocessed_reads.sample.tsv"); do
	dir=$(dirname $tsv)
	sed "s|/nfs3_ib/nfs-ip34||g" ${tsv} > ${tsv}.fast
	sed -i "s|/net/nfs-ip34||g" ${tsv}.fast
	sed -i "s|${dir}|/fast2/def-ilafores/preproc|g" ${tsv}.fast
done

ml apptainer
k2_std=/dev/shm/k2_standard_20241228 
k2_gtdb=/dev/shm/k2_gtdb_genome_reps_20241109
k2_local=$MC/scripts/kraken_local.sh

# RUN KRAKEN ON /fast2/
bash $k2_local --tsv ${TSV}.fast --confidence 0.10 --output $MC/$DATASET/KB10 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${TSV}.fast --confidence 0.45 --output $MC/$DATASET/KB45 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${TSV}.fast --confidence 0.90 --output $MC/$DATASET/KB90 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${TSV}.fast --confidence 0.10 --output $MC/$DATASET/KB10_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${TSV}.fast --confidence 0.45 --output $MC/$DATASET/KB45_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${TSV}.fast --confidence 0.90 --output $MC/$DATASET/KB90_GTDB --kraken_db $k2_gtdb --threads 24

# Record % sequence classification
cd "$MC" && find . -name '*.kreport' ! -name '*bracken*' -exec awk '
FNR==1 {unclassified=$2; next} 
FNR==2 {classified=$2; rate=classified/(unclassified+classified); 
         printf "%s\t%d\t%d\t%.5f\n", FILENAME, classified, unclassified, rate}
' {} + > kraken_classification_rate.tsv

# Check completion status
check_output 'KB10 KB45 KB90 KB10_GTDB KB45_GTDB KB90_GTDB' $DATASET _bracken_S.MPA.TXT

database="KB90"
missing_KB=$(grep -n -v -f <(ls $DATASET/$database/*/*/*_bracken_S.MPA.TXT | awk -F'/' '{print $3}' | sed 's/_profile\.txt//') $DATASET/preproc/preprocessed_reads.sample.tsv | cut -f1 -d: | tr '\n' ','); echo $missing_KB
sbatch --array="$missing_KB" $MC/$DATASET/$database/taxonomic_profile.samples.slurm.sh $DATASET "\$${dataset}_TSV"

# Once completely done, remove heavy files from kraken out 
rm */KB*/*/*_taxonomy_nt */KB*/*/*/*.bracken */KB*/*/*bugs_list.MPA.TXT */KB*/*/*/*_temp.MPA.TXT */KB*/*/*/*_bracken_[^S].MPA.TXT -r */KB*/*/*_kronagrams -r */*/.throttle/

################
# MetaPhlAn4 ###
################

metaphlan="bash $ANCHOR/$ILL_PIPELINES/generateslurm_taxonomic_abundance.metaphlan.sh \
	--slurm_log $ANCHOR/$MC/logs --slurm_walltime 24:00:00 --slurm_threads 24 --slurm_mem 30G"

# Generate SLURM scripts https://github.com/jflucier/ILL_pipelines/blob/main/generateslurm_taxonomic_abundance.metaphlan.sh
# 2022 database

$metaphlan --sample_tsv $ANCHOR/$TSV --db $FAST/metaphlan3_db/mpa_v30_CHOCOPhlAn_201901 --out $ANCHOR/$MC/$DATASET/MPA_db2019
sbatch --array=1-"$N_SAMPLES" $ANCHOR/$MC/$DATASET/MPA_db2019/metaphlan.slurm.sh

$metaphlan --sample_tsv $ANCHOR/$TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $ANCHOR/$MC/$DATASET/MPA_db2022
sbatch --array=1-"$N_SAMPLES" $ANCHOR/$MC/$DATASET/MPA_db2022/metaphlan.slurm.sh

$metaphlan --sample_tsv $ANCHOR/$TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $ANCHOR/$MC/$DATASET/MPA_db2023
sbatch --array=1-"$N_SAMPLES" $ANCHOR/$MC/$DATASET/MPA_db2023/metaphlan.slurm.sh

# Rerun missing MPA
database="MPA_db2023"
missing_MPA=$(grep -n -v -f <(ls $DATASET/$database/*/*_profile.txt | awk -F'/' '{print $3}' | sed 's/_profile\.txt//') $DATASET/preproc/preprocessed_reads.sample.tsv | cut -f1 -d: | tr '\n' ','); echo $missing_MPA
sbatch --array="$missing_MPA" $MC/$DATASET/$database/metaphlan.slurm.sh $DATASET "\$${dataset}_TSV"

# Check completion status
check_output 'MPA_db2022 MPA_db2023' $DATASET _profile.txt

# Remove bowtie indexes
rm */MPA_db*/*/*.bowtie2.txt
rm */*/.throttle -r

#####################
# Sourmash gather ###
#####################

sbatch --mem=120G -n 24 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "genbank-2022.03"
sbatch --mem=31G -n 24 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "gtdb-rs214-rep"
sbatch --mem=80G -n 16 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "gtdb-rs214-full"
sbatch --mem=31G -n 24 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "gtdb-rs220-rep"
sbatch --mem=60G -n 12 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "/fast2/def-ilafores/refseq_genomes/refseq-229.k31.sig"

# Check completion status
check_output 'gtdb-rs214-rep gtdb-rs214-full genbank-2022.03 gtdb-reps-rs220' $DATASET _gather.csv

# Extract Sourmash lineage subset
for i in P19_Saliva P19_Gut PD AD_Skin Moss RA_Gut Bee Olive NAFLD; do
for j in SM_gtdb-rs214-full SM_gtdb-rs214-rep SM_genbank-2022.03 SM_gtdb-rs220-rep; do
	db=$(echo "$j" | cut -d'_' -f2 | cut -d'-' -f1,2)
	ver=$(echo "$j" | cut -d'_' -f3)

	if [[ ! -f $i/$j/${j}_lineages.csv ]]; then
        cat $i/$j/*_gather.csv | cut -d, -f10 | tail -n+2 | awk '{print $1}' | sed 's/"//; s/^[^_]*_//; s/\.[^.]*$//' | grep -v name | sort -u | grep -Fhf - $ILAFORES/ref_dbs/sourmash_db/${db}*${ver}*.lineages.csv > $i/$j/${j}_lineages.csv
	fi
done
done

## surplus taxa in sourmash rs220 index
cat $MC/P19_Gut/Sourmash/*rs220*_gather.csv | cut -d, -f10 | tail -n+2 | \
	awk '{print $1}' | sed 's/"//' | sort -u > found_taxa.tsv # 8070 taxa

# of which 1191 are not in the species reps lineage file
grep -v -f <(cut -f1 $ILAFORES/ref_dbs/sourmash_db/bac120_taxonomy_r220.tsv | sed 's/^[^_]*_//' | sort -u) found_taxa.tsv | wc






## Remove and Redownload corrupted samples:
remove_these=($(sed -n "$(echo $missing_samples | sed 's/,/p;/g' | sed 's/;$//')" $DATASET/preproc/preprocessed_reads.sample.tsv | awk '{print $1}' ))
for sample in "${remove_these[@]}"; do
	rm PD/raw/${sample}_?.fastq.gz
done
# download!
cd PD/raw
grep -f <(printf "%s\n" "${remove_these[@]}") "../ena-file-download-read_run-PRJNA834801-fastq_ftp-20250310-1443.sh" | bash
cd ../../






# Find missing line numbers
grep -n -v -f <(ls NAFLD/SM_genbank_202203/ | sed 's/_.*//' | sort | uniq) NAFLD/preproc/preprocessed_reads.sample.tsv | awk -F: '{print $1}' | paste -sd,

:>$MC/NAFLD/raw/samples_to_process.tsv
while read -r p; do 
SRR=$(echo $p | awk '{print $1}')
f1=$(find NAFLD/raw -type f -name "${SRR}_1.fastq*" -exec realpath {} \;)
f2=$(find NAFLD/raw -name "${SRR}_2.fastq*" -exec realpath {} \;)
SAM=$(echo $p | awk '{print $14}')
if [ ! -z "$f1" ]; then 
	echo -e "$SAM\t$f1\t$f2" >> $MC/NAFLD/raw/samples_to_process.tsv
fi; done < NAFLD/raw/ENA_report.tsv 


:> $MC/NAFLD/raw/samples_to_process.tsv
while read p; do 
SAM=$(echo $p | cut -f1 -d' ')
if [[ "$SAM" == "run_accession" ]]; then 
continue
fi
f1=$(find $MC/NAFLD/raw -name "${SAM}_1.fastq.gz")
f2=$(find $MC/NAFLD/raw -name "${SAM}_2.fastq.gz")
echo -e "$SAM\t${f1}\t${f2}" >> $MC/NAFLD/raw/samples_to_process.tsv
done < $MC/NAFLD/raw/metadata.tsv 


### Update kraken db
bracken="bash $ILL_PIPELINES/generateslurm_taxonomic_profile.sample.sh \
	--kraken_db $ILAFORES/ref_dbs/kraken2_dbs/kraken2_PlusPFP_202202 \
	--slurm_log $MC/logs --slurm_walltime 72:00:00 --slurm_threads 48 --slurm_mem 250G"

