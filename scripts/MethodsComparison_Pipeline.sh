export MC=$ILAFORES/analysis/MethodsComparison
#export ILL_PIPELINES=$ILAFORES/analysis/MethodsComparison/ILL_pipelines
cd $MC
source scripts/myFunctions.sh

# Create <dataset>_TSV and NUM_<dataset> variables
# Choose the one to work with :
dataset_variables "P19_Saliva" "$PR19/Saliva/preproc/preprocessed_reads.sample.tsv"
dataset_variables "P19_Gut" "$PR19/Feces/preproc/preprocessed_reads.sample.tsv"
dataset_variables "Moss" "$MOSS/preproc/preprocessed_reads.sample.tsv"
dataset_variables "NAFLD" "$MC/data/NAFLD/preproc/preprocessed_reads.sample.tsv"
dataset_variables "AD_Skin" "$MC/data/AD_Skin/preproc/preprocessed_reads.sample.tsv"
dataset_variables "PD" "$MC/data/PD/preproc/preprocessed_reads.sample.tsv"
dataset_variables "Bee" "$MC/data/Bee/preproc/preprocessed_reads.sample.tsv"
dataset_variables "Olive" "$MC/data/Olive/preproc/preprocessed_reads.sample.tsv"
dataset_variables "RA_Gut" "$MC/data/RA_Gut/preproc/preprocessed_reads.sample.tsv"

######################
# QC #################
######################

mkdir -p  $DATASET_PATH/raw

# Download from ENA tsv file with grep subset
samples=($(grep "WGS"  $DATASET_PATH/raw/filereport_*_tsv.txt | cut -f1  | grep -v 'run'))
cd /fast2/def-ilafores/Olive/raw
grep -E "$(IFS="|"; echo "${samples[*]}")" " $DATASET_PATH/raw/" | bash

# Create TSV
:>  $DATASET_PATH/raw/samples_to_process.tsv
for sample in "${samples[@]}"; do
	fq1=$(find $ANCHOR/fast2/def-ilafores/$DATASET/raw -type f -name "${sample}_1.fastq.gz")
	fq2=$(find $ANCHOR/fast2/def-ilafores/$DATASET/raw -type f -name "${sample}_2.fastq.gz")
	echo -e "${sample}\t${fq1}\t${fq2}" >>  $DATASET_PATH/raw/samples_to_process.tsv
done

#Process samples
mkdir -p /fast2/def-ilafores/$DATASET/preproc/logs
bash $ANCHOR/$ILAFORES/programs/ILL_pipelines/generateslurm_preprocess.kneaddata.sh \
	--sample_tsv $ANCHOR/ $DATASET_PATH/raw/samples_to_process.tsv \
	--out $ANCHOR/fast2/def-ilafores/$DATASET/preproc \
	--trimmomatic_options "SLIDINGWINDOW:4:20 MINLEN:50" \
	--db $FAST/host_genomes/GRCh38_index/grch38_1kgmaj \
	--slurm_mem 30G --slurm_threads 24
# correct script as the shell command needs to be anchored!
cp -r /fast2/def-ilafores/$DATASET/preproc/*  $DATASET_PATH/preproc/

# Find arrays of missing samples:
missing_samples=$(grep -n -v -f <(ls  $DATASET_PATH/preproc/*/*_1.fastq.gz | awk -F'/' '{print $3}')  $DATASET_PATH/raw/samples_to_process.tsv | cut -f1 -d: | tr '\n' ','); echo $missing_samples

rm -r  $DATASET_PATH/preproc/.throttle
sbatch --array="$missing_samples" /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/PD/preproc/preprocess.kneaddata.slurm.sh

# Sequence counts by sample
ml seqkit 
:>$MC/Out/stats/all_read_counts_raw.tsv

fastq_files=($(find ./Bee/preproc -maxdepth 2 -type f -name "*paired_1.fastq*"))
total_files=${#fastq_files[@]}

start=0
chunk_size=40

while [ $start -lt $total_files ]; do
    seqkit stats ${fastq_files[@]:$start:$chunk_size} --threads $chunk_size --skip-err >> $MC/Out/stats/all_read_counts_raw.tsv
    start=$((start + chunk_size))
done

grep DNA $MC/Out/stats/all_read_counts_raw.tsv | sponge $MC/Out/stats/all_read_counts_raw.tsv

# Replace first column with dataset and sample IDs :
awk '{
    path = $1;
    rest = substr($0, length(path) + 2);
    split(path, parts, "/");
    group = parts[2];
    sample = parts[4];
    print group " " sample " " rest;
}' OFS='\t' $MC/Out/stats/all_read_counts_raw.tsv > $MC/Out/stats/all_read_counts.tsv

# Summarise mean ± SD by dataset in awk lol just do R gobless deepseek
awk -F' ' '
{
    gsub(/,/, "", $5);  # Remove commas from the fifth column
    group = $1;
    value = $5 + 0;
    
    # Initialize min/max for new groups
    if (!(group in count)) {
        min[group] = value;
        max[group] = value;
    }
    
    # Update min/max
    if (value < min[group]) min[group] = value;
    if (value > max[group]) max[group] = value;
    
    # Track sum and sumsq for mean/stddev
    count[group]++;
    sum[group] += value;
    sumsq[group] += value * value;
}
END {
    for (g in count) {
        mean = sum[g] / count[g];
        stddev = sqrt((sumsq[g] - sum[g]^2 / count[g]) / count[g]);
        printf "%s\t%.2f ± %.2f [%d, %d]\n", g, mean, stddev, min[g], max[g];
    }
}' $MC/Out/stats/all_read_counts.tsv > $MC/Out/stats/sample_count_summary.tsv

################
# mOTUs ########
################

# Custom SLURM script
sbatch --array=1-"$N_SAMPLES" $MC/scripts/motus_SLURM.sh  $DATASET_PATH $TSV.fast

# Check completion status
check_output 'MOTUS'  $DATASET_PATH _profile.txt
 
# Rerun missing MOTUS
rm  $DATASET_PATH/MOTUS/_profile.txt # not sure why that appears
missing_motus=$(grep -n -v -f <(ls " $DATASET_PATH/MOTUS/"*_profile.txt | awk -F'/' '{print $3}' | sed 's/_profile\.txt//') "$(eval echo \$${dataset}_TSV)" | cut -f1 -d: | tr '\n' ','); echo "$missing_motus"
sbatch --array="$missing_motus" $MC/scripts/motus_SLURM.sh  $DATASET_PATH "$(eval echo \$${dataset}_TSV)"

####################
# Kraken/bracken ###
# on /fast2/ local #
####################

# Copy fastqs to /fast2
for dir in $(find $MC/data -maxdepth 3 -type d -name 'preproc'); do 
nice -n10 ionice -c2 -n7 rclone copy $dir /fast2/def-ilafores/preproc --transfers 16 --checkers 4 --modify-window 5s --fast-list --no-update-modtime --retries 3 --retries-sleep 30s --low-level-retries 1 -v -L --size-only --exclude "*contam*";
done

# rearrange tsvs to point to new path
for tsv in $(find $ILAFORES/analysis/ -maxdepth 5 -name "preprocessed_reads.sample.tsv"); do
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
bash $k2_local --tsv ${TSV}.fast --confidence 0.10 --output  $DATASET_PATH/KB10 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${TSV}.fast --confidence 0.45 --output  $DATASET_PATH/KB45 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${TSV}.fast --confidence 0.90 --output  $DATASET_PATH/KB90 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${TSV}.fast --confidence 0.10 --output  $DATASET_PATH/KB10_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${TSV}.fast --confidence 0.45 --output  $DATASET_PATH/KB45_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${TSV}.fast --confidence 0.90 --output  $DATASET_PATH/KB90_GTDB --kraken_db $k2_gtdb --threads 24

# Record % sequence classification
mkdir -p Out/classification_rates
cd "$MC" && find . -name '*.kreport' ! -name '*bracken*' -exec awk '
FNR==1 {unclassified=$2; next} 
FNR==2 {classified=$2; rate=classified/(unclassified+classified); 
         printf "%s\t%d\t%d\t%.5f\n", FILENAME, classified, unclassified, rate}
' {} + > Out/classification_rates/kraken_classification_rate.tsv

# Check completion status
check_output 'KB10 KB45 KB90 KB10_GTDB KB45_GTDB KB90_GTDB' "$DATASET_PATH" _bracken_S.MPA.TXT

database="KB90"
missing_KB=$(grep -n -v -f <(ls  $DATASET_PATH/$database/*/*/*_bracken_S.MPA.TXT | awk -F'/' '{print $3}' | sed 's/_profile\.txt//')  $DATASET_PATH/preproc/preprocessed_reads.sample.tsv | cut -f1 -d: | tr '\n' ','); echo $missing_KB
sbatch --array="$missing_KB"  $DATASET_PATH/$database/taxonomic_profile.samples.slurm.sh  $DATASET_PATH "\$${dataset}_TSV"

# Once completely done, remove heavy files from kraken out 
rm */KB*/*/*_taxonomy_nt */KB*/*/*/*.bracken */KB*/*/*bugs_list.MPA.TXT */KB*/*/*/*_temp.MPA.TXT */KB*/*/*/*_bracken_[^S].MPA.TXT -r */KB*/*/*_kronagrams -r */*/.throttle/

################
# MetaPhlAn4 ###
################

metaphlan="bash $ANCHOR/$ILL_PIPELINES/generateslurm_taxonomic_abundance.metaphlan.sh \
	--slurm_log $ANCHOR/$MC/logs --slurm_walltime 24:00:00 --slurm_threads 24 --slurm_mem 30G"

# Generate SLURM scripts https://github.com/jflucier/ILL_pipelines/blob/main/generateslurm_taxonomic_abundance.metaphlan.sh
# 2022 database

$metaphlan --sample_tsv $ANCHOR/$TSV --db $FAST/metaphlan3_db/mpa_v30_CHOCOPhlAn_201901 --out $ANCHOR/ $DATASET_PATH/MPA_db2019
sbatch --array=1-"$N_SAMPLES" $ANCHOR/ $DATASET_PATH/MPA_db2019/metaphlan.slurm.sh

$metaphlan --sample_tsv $ANCHOR/$TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $ANCHOR/ $DATASET_PATH/MPA_db2022
sbatch --array=1-"$N_SAMPLES" $ANCHOR/ $DATASET_PATH/MPA_db2022/metaphlan.slurm.sh

$metaphlan --sample_tsv $ANCHOR/$TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $ANCHOR/ $DATASET_PATH/MPA_db2023
sbatch --array=1-"$N_SAMPLES" $ANCHOR/ $DATASET_PATH/MPA_db2023/metaphlan.slurm.sh

# Rerun missing MPA
database="MPA_db2023"
missing_MPA=$(grep -n -v -f <(ls  $DATASET_PATH/$database/*/*_profile.txt | awk -F'/' '{print $3}' | sed 's/_profile\.txt//')  $DATASET_PATH/preproc/preprocessed_reads.sample.tsv | cut -f1 -d: | tr '\n' ','); echo $missing_MPA
sbatch --array="$missing_MPA"  $DATASET_PATH/$database/metaphlan.slurm.sh  $DATASET_PATH "\$${dataset}_TSV"

# Check completion status
check_output 'MPA_db2022 MPA_db2023'  $DATASET_PATH _profile.txt

# Remove bowtie indexes
rm */MPA_db*/*/*.bowtie2.txt
rm */*/.throttle -r

#####################
# Sourmash gather ###
#####################

sbatch --mem=120G -n 24 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "genbank-2022.03"
sbatch --mem=31G -n 24 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "gtdb-rs214-rep"
sbatch --mem=80G -n 16 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "gtdb-rs214-full"
sbatch --mem=31G -n 24 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "/home/def-ilafores/ref_dbs/sourmash_db/gtdb-rs220-reps-k31.zip"
sbatch --mem=120G -n 24 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "/fast2/def-ilafores/refseq_genomes/all_sig_refseq-09-April-2025.zip"
sbatch --mem=31G -n 24 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "/home/def-ilafores/analysis/boreal_moss/genome_sketches/gtdb-rs214-rep-MAGs.sbt.zip"
sbatch --mem=31G -n 1 --array=1-"$N_SAMPLES" $MC/scripts/gather_SLURM_fast.sh "$DATASET" "$TSV.fast" "$ILAFORES/ref_dbs/sourmash_db/RefSeq_20250528.k31.sbt.zip" # Here using 1 see gather_SLURM_fast.sh script; pending sourmash 4.9 update to remove 

# Check completion status
check_output 'RefSeq_20250528'  $DATASET_PATH _gather.csv

# Extract Sourmash lineage subset
for i in P19_Saliva P19_Gut PD AD_Skin Moss RA_Gut Bee Olive NAFLD; do
for j in $(find data/$i -maxdepth 1 -type d -name 'SM_*'); do
	db_name=$(basename $j)
	db=$(echo "$db_name" | cut -d'_' -f2 | cut -d'-' -f1,2)
	ver=$(echo "$db_name" | cut -d'_' -f3)

	if [[ ! -f $j/${db_name}_lineages.csv ]]; then
        cat $j/*_gather.csv | cut -d, -f10 | tail -n+2 | awk '{print $1}' | sed 's/"//; s/^[^_]*_//; s/\.[^.]*$//' | grep -v name | sort -u | grep -Fhf - $ILAFORES/ref_dbs/sourmash_db/${db}*${ver}*.lineages.csv > $j/${db_name}_lineages.csv
	fi
done
done

# Compute classification rates from sourmash gather output : 
mkdir -p Out/classification_rates
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

## surplus taxa in sourmash rs220 index
mkdir tmp
cat $MC/data/P19_Gut/Sourmash/*rs220*_gather.csv | cut -d, -f10 | tail -n+2 | \
	awk '{print $1}' | sed 's/"//' | sort -u > tmp/found_taxa.tsv # 8070 taxa

# of which 1191 are not in the species reps lineage file
grep -v -f <(cut -f1 $ILAFORES/ref_dbs/sourmash_db/bac120_taxonomy_r220.tsv | sed 's/^[^_]*_//' | sort -u) tmp/found_taxa.tsv | wc


## Remove and Redownload corrupted samples:
remove_these=($(sed -n "$(echo $missing_samples | sed 's/,/p;/g' | sed 's/;$//')"  $DATASET_PATH/preproc/preprocessed_reads.sample.tsv | awk '{print $1}' ))
for sample in "${remove_these[@]}"; do
	rm data/PD/raw/${sample}_?.fastq.gz
done
# download!
cd data/PD/raw
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

