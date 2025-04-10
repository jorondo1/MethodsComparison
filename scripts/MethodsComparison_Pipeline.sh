export MC=$ILAFORES/analysis/MethodsComparison
export ILL_PIPELINES=$ILAFORES/analysis/MethodsComparison/ILL_pipelines
cd $MC
source scripts/myFunctions.sh
export SALIVA_TSV=$PR19/Saliva/preproc/preprocessed_reads.sample.tsv
export FECES_TSV=$PR19/Feces/preproc/preprocessed_reads.sample.tsv
export MOSS_TSV=$MOSS/preproc/preprocessed_reads.sample.tsv
export NAFLD_TSV=$MC/NAFLD/preproc/preprocessed_reads.sample.tsv
export AD_Skin_TSV=$MC/AD_Skin/preproc/preprocessed_reads.sample.tsv
export PD_TSV=$MC/PD/preproc/preprocessed_reads.sample.tsv
export NUM_P19_Saliva=$(wc $SALIVA_TSV | awk '{print $1}')
export NUM_P19_Gut=$(wc $FECES_TSV | awk '{print $1}')
export NUM_Moss=$(wc $MOSS_TSV | awk '{print $1}')
export NUM_NAFLD=$(wc $NAFLD_TSV | awk '{print $1}')
export NUM_AD_Skin=$(wc $AD_Skin_TSV | awk '{print $1}')
export NUM_PD=$(wc $PD_TSV | awk '{print $1}')
export DATASETS="P19_Saliva P19_Gut Moss NAFLD AD_Skin PD"

######################
# QC #################
######################

#Process fecal samples
dataset="Bee"
mkdir -p $dataset/preproc
bash $ILAFORES/programs/ILL_pipelines/generateslurm_preprocess.kneaddata.sh \
	--sample_tsv $MC/$dataset/raw/samples_to_process.tsv \
	--out $MC/$dataset/preproc \
	--trimmomatic_options "SLIDINGWINDOW:4:20 MINLEN:50" \
	--db $FAST/host_genomes/GRCh38_index/grch38_1kgmaj \
	--slurm_mem 120G --slurm_threads 24
# correct script as the shell command needs to be anchored!

# Find arrays of missing samples:
missing_samples=$(grep -n -v -f <(ls $dataset/preproc/*/*_1.fastq.gz | awk -F'/' '{print $3}') $dataset/raw/samples_to_process.tsv | cut -f1 -d: | tr '\n' ','); echo $missing_samples

rm -r $dataset/preproc/.throttle
sbatch --array="$missing_samples" /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/PD/preproc/preprocess.kneaddata.slurm.sh

############
# mOTUs ####
# Custom SLURM script 
sbatch --array=1-"$NUM_P19_Saliva" $MC/scripts/motus_SLURM.sh P19_Saliva $SALIVA_TSV
sbatch --array=1-"$NUM_P19_Gut" $MC/scripts/motus_SLURM.sh P19_Gut $FECES_TSV
sbatch --array=1-"$NUM_Moss" $MC/scripts/motus_SLURM.sh Moss $MOSS_TSV
sbatch --array=1-"$NUM_NAFLD" $MC/scripts/motus_SLURM.sh NAFLD $NAFLD_TSV
sbatch --array=1-"$NUM_AD_Skin" $MC/scripts/motus_SLURM.sh AD_Skin $AD_Skin_TSV
sbatch --array=1-"$NUM_PD" $MC/scripts/motus_SLURM.sh PD $PD_TSV

# Check completion status
check_output 'MOTUS' 'PD' _profile.txt

# Rerun missing MOTUS
dataset="PD"
rm $dataset/MOTUS/_profile.txt # not sure why that appears
missing_motus=$(grep -n -v -f <(ls "$dataset/MOTUS/"*_profile.txt | awk -F'/' '{print $3}' | sed 's/_profile\.txt//') "$(eval echo \$${dataset}_TSV)" | cut -f1 -d: | tr '\n' ','); echo "$missing_motus"
sbatch --array="$missing_motus" $MC/scripts/motus_SLURM.sh $dataset "$(eval echo \$${dataset}_TSV)"

####################
# Kraken/bracken ###

## Run kraken local ip34 using /fast2
# Copy fastqs to /fast2

for dir in $(find $MC -maxdepth 3 -type d -name 'preproc'); do 
nice -n10 ionice -c2 -n7 rclone copy $dir /fast2/def-ilafores/preproc --transfers 16 --checkers 4 --modify-window 5s --fast-list --no-update-modtime --retries 3 --retries-sleep 30s --low-level-retries 1 -v -L --size-only --exclude "*contam*";
done

# rearrange tsvs to point to new path
for tsv in $(find $ILAFORES/analysis/ -name "preprocessed_reads.sample.tsv"); do
	dir=$(dirname $tsv)
	sed "s|/nfs3_ib/nfs-ip34||g" ${tsv} > ${tsv}.fast
	sed -i "s|/net/nfs-ip34||g" ${tsv}.fast
	sed -i "s|${dir}|/fast2/def-ilafores/preproc|g" ${tsv}.fast
done

ml apptainer
k2_std=/dev/shm/k2_standard_20241228 
k2_gtdb=/dev/shm/k2_gtdb_genome_reps_20241109
k2_local=$MC/scripts/kraken_local.sh

# PD
bash $k2_local --tsv ${PD_TSV}.fast --confidence 0.10 --output $MC/PD/KB10 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${PD_TSV}.fast --confidence 0.45 --output $MC/PD/KB45 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${PD_TSV}.fast --confidence 0.90 --output $MC/PD/KB90 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${PD_TSV}.fast --confidence 0.10 --output $MC/PD/KB10_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${PD_TSV}.fast --confidence 0.45 --output $MC/PD/KB45_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${PD_TSV}.fast --confidence 0.90 --output $MC/PD/KB90_GTDB --kraken_db $k2_gtdb --threads 24

# Moss
bash $k2_local --tsv ${MOSS_TSV}.fast --confidence 0.10 --output $MC/Moss/KB10 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${MOSS_TSV}.fast --confidence 0.45 --output $MC/Moss/KB45 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${MOSS_TSV}.fast --confidence 0.90 --output $MC/Moss/KB90 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${MOSS_TSV}.fast --confidence 0.10 --output $MC/Moss/KB10_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${MOSS_TSV}.fast --confidence 0.45 --output $MC/Moss/KB45_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${MOSS_TSV}.fast --confidence 0.90 --output $MC/Moss/KB90_GTDB --kraken_db $k2_gtdb --threads 24

# NAFLD
bash $k2_local --tsv ${NAFLD_TSV}.fast --confidence 0.10 --output $MC/NAFLD/KB10 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${NAFLD_TSV}.fast --confidence 0.45 --output $MC/NAFLD/KB45 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${NAFLD_TSV}.fast --confidence 0.90 --output $MC/NAFLD/KB90 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${NAFLD_TSV}.fast --confidence 0.10 --output $MC/NAFLD/KB10_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${NAFLD_TSV}.fast --confidence 0.45 --output $MC/NAFLD/KB45_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${NAFLD_TSV}.fast --confidence 0.90 --output $MC/NAFLD/KB90_GTDB --kraken_db $k2_gtdb --threads 24

#AD_SKIN
bash $k2_local --tsv ${AD_Skin_TSV}.fast --confidence 0.10 --output $MC/AD_Skin/KB10 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${AD_Skin_TSV}.fast --confidence 0.45 --output $MC/AD_Skin/KB45 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${AD_Skin_TSV}.fast --confidence 0.90 --output $MC/AD_Skin/KB90 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${AD_Skin_TSV}.fast --confidence 0.10 --output $MC/AD_Skin/KB10_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${AD_Skin_TSV}.fast --confidence 0.45 --output $MC/AD_Skin/KB45_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${AD_Skin_TSV}.fast --confidence 0.90 --output $MC/AD_Skin/KB90_GTDB --kraken_db $k2_gtdb --threads 24

# P19_Feces
bash $k2_local --tsv ${FECES_TSV}.fast --confidence 0.10 --output $MC/P19_Feces/KB10 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${FECES_TSV}.fast --confidence 0.45 --output $MC/P19_Feces/KB45 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${FECES_TSV}.fast --confidence 0.90 --output $MC/P19_Feces/KB90 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${FECES_TSV}.fast --confidence 0.10 --output $MC/P19_Feces/KB10_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${FECES_TSV}.fast --confidence 0.45 --output $MC/P19_Feces/KB45_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${FECES_TSV}.fast --confidence 0.90 --output $MC/P19_Feces/KB90_GTDB --kraken_db $k2_gtdb --threads 24

# P19_Saliva
bash $k2_local --tsv ${SALIVA_TSV}.fast --confidence 0.10 --output $MC/P19_Saliva/KB10 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${SALIVA_TSV}.fast --confidence 0.45 --output $MC/P19_Saliva/KB45 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${SALIVA_TSV}.fast --confidence 0.90 --output $MC/P19_Saliva/KB90 --kraken_db $k2_std --threads 24
bash $k2_local --tsv ${SALIVA_TSV}.fast --confidence 0.10 --output $MC/P19_Saliva/KB10_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${SALIVA_TSV}.fast --confidence 0.45 --output $MC/P19_Saliva/KB45_GTDB --kraken_db $k2_gtdb --threads 24
bash $k2_local --tsv ${SALIVA_TSV}.fast --confidence 0.90 --output $MC/P19_Saliva/KB90_GTDB --kraken_db $k2_gtdb --threads 24

# Record % sequence classification
cd $MC; find . -name '*.kreport' ! -name '*bracken*' -exec awk '
FNR==1 {unclassified=$2; next} 
FNR==2 {classified=$2; printf "%s\t%.5f\n", FILENAME, (classified/(unclassified+classified))}
' {} + > kraken_classification_rate.tsv

# Check completion status
check_output 'KB10 KB45 KB90' 'PD' _bracken_S.MPA.TXT

dataset="PD"
database="KB90"
missing_KB=$(grep -n -v -f <(ls $dataset/$database/*/*/*_bracken_S.MPA.TXT | awk -F'/' '{print $3}' | sed 's/_profile\.txt//') $dataset/preproc/preprocessed_reads.sample.tsv | cut -f1 -d: | tr '\n' ','); echo $missing_KB
sbatch --array="$missing_KB" $MC/$dataset/$database/taxonomic_profile.samples.slurm.sh $dataset "\$${dataset}_TSV"

# Once completely done, remove taxonomy_nt files from kraken out 
rm */KB*/*/*_taxonomy_nt # heavy
rm */KB*/*/*/*.bracken # useless format
rm */KB*/*/*/*.kreport # we want the bracken out, not kraken
rm */KB*/*/*bugs_list.MPA.TXT # no use
rm */KB*/*/*/*_temp.MPA.TXT # temp files
rm */KB*/*/*/*_bracken_[^S].MPA.TXT # Other levels of classification
rm -r */KB*/*/*_kronagrams # no use
rm -r */*/.throttle/

################
# MetaPhlAn4 ###
metaphlan="bash $ILL_PIPELINES/generateslurm_taxonomic_abundance.metaphlan.sh \
	--slurm_log $MC/logs --slurm_walltime 24:00:00 --slurm_threads 24 --slurm_mem 30G"

# Generate SLURM scripts https://github.com/jflucier/ILL_pipelines/blob/main/generateslurm_taxonomic_abundance.metaphlan.sh
# 2022 database
$metaphlan --sample_tsv $SALIVA_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/P19_Saliva/MPA_db2022
$metaphlan --sample_tsv $SALIVA_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/P19_Saliva/MPA_db2023
sbatch --array=1-"$NUM_P19_Saliva" $MC/P19_Saliva/MPA_db2022/metaphlan.slurm.sh
sbatch --array=1-"$NUM_P19_Saliva" $MC/P19_Saliva/MPA_db2023/metaphlan.slurm.sh

$metaphlan --sample_tsv $FECES_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/P19_Gut/MPA_db2022
$metaphlan --sample_tsv $FECES_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/P19_Gut/MPA_db2023
sbatch --array=1-"$NUM_P19_Gut" $MC/P19_Gut/MPA_db2022/metaphlan.slurm.sh
sbatch --array=1-"$NUM_P19_Gut" $MC/P19_Gut/MPA_db2023/metaphlan.slurm.sh

$metaphlan --sample_tsv $MOSS_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/Moss/MPA_db2022
$metaphlan --sample_tsv $MOSS_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/Moss/MPA_db2023
sbatch --array=1-"$NUM_Moss" $MC/Moss/MPA_db2022/mfetaphlan.slurm.sh
sbatch --array=1-"$NUM_Moss" $MC/Moss/MPA_db2023/metaphlan.slurm.sh

$metaphlan --sample_tsv $NAFLD_TSV --db $FAST/metaphlan3_db/mpa_v30_CHOCOPhlAn_201901 --out $MC/NAFLD/MPA_db2019
$metaphlan --sample_tsv $NAFLD_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/NAFLD/MPA_db2022
$metaphlan --sample_tsv $NAFLD_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/NAFLD/MPA_db2023
sbatch --array=1-"$NUM_NAFLD" $MC/NAFLD/MPA_db2019/metaphlan.slurm.sh
sbatch --array=1-"$NUM_NAFLD" $MC/NAFLD/MPA_db2022/metaphlan.slurm.sh
sbatch --array=1-"$NUM_NAFLD" $MC/NAFLD/MPA_db2023/metaphlan.slurm.sh

$metaphlan --sample_tsv $AD_Skin_TSV --db $FAST/metaphlan3_db/mpa_v30_CHOCOPhlAn_201901 --out $MC/AD_Skin/MPA_db2019
$metaphlan --sample_tsv $AD_Skin_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/AD_Skin/MPA_db2022
$metaphlan --sample_tsv $AD_Skin_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/AD_Skin/MPA_db2023
sbatch --array=1-"$NUM_AD_Skin" $MC/AD_Skin/MPA_db2019/metaphlan.slurm.sh
sbatch --array=1-"$NUM_AD_Skin" $MC/AD_Skin/MPA_db2022/metaphlan.slurm.sh
sbatch --array=1-"$NUM_AD_Skin" $MC/AD_Skin/MPA_db2023/metaphlan.slurm.sh

$metaphlan --sample_tsv $PD_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/PD/MPA_db2022
$metaphlan --sample_tsv $PD_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/PD/MPA_db2023
sbatch --array=1-"$NUM_PD" $MC/PD/MPA_db2022/metaphlan.slurm.sh
sbatch --array=1-"$NUM_PD" $MC/PD/MPA_db2023/metaphlan.slurm.sh

#singularity exec -e -B $ILAFORES:$ILAFORES -B $FAST:$FAST $ILAFORES/programs/ILL_pipelines/containers/humann.3.6.sif metaphlan \
#	-t rel_ab --input_type fastq --unclassified_estimation --mpa3 \
#	--bowtie2db $FAST/metaphlan3_db \
#	-x 'mpa_v30_CHOCOPhlAn_201901' \
#	--nproc 72 $MC/AD_Skin/preproc/SRR20002471/SRR20002471_1.fastq.gz $PWD/test_profile.txt 
#	

# Rerun missing MPA
dataset="PD"
database="MPA_db2023"
missing_MPA=$(grep -n -v -f <(ls $dataset/$database/*/*_profile.txt | awk -F'/' '{print $3}' | sed 's/_profile\.txt//') $dataset/preproc/preprocessed_reads.sample.tsv | cut -f1 -d: | tr '\n' ','); echo $missing_MPA
sbatch --array="$missing_MPA" $MC/$dataset/$database/metaphlan.slurm.sh $dataset "\$${dataset}_TSV"

# Check completion status
check_output 'MPA_db2022 MPA_db2023 MPA_db2019' 'PD' _profile.txt

# Remove bowtie indexes
rm */MPA_db*/*/*.bowtie2.txt
rm */*/.throttle -r

#####################
# Sourmash gather ###

sbatch --mem=60G -n 16 --array=1-"$NUM_P19_Saliva" $MC/scripts/gather_SLURM.sh "P19_Saliva" $SALIVA_TSV "genbank-2022.03"
sbatch --mem=31G -n 24 --array=1-"$NUM_P19_Saliva" $MC/scripts/gather_SLURM.sh "P19_Saliva" $SALIVA_TSV "gtdb-rs214-rep"
sbatch --mem=80G -n 16 --array=1-"$NUM_P19_Saliva" $MC/scripts/gather_SLURM.sh "P19_Saliva" $SALIVA_TSV "gtdb-rs214-full"

sbatch --mem=80G -n 16 --array=1-"$NUM_P19_Gut" $MC/scripts/gather_SLURM.sh "P19_Gut" $FECES_TSV "genbank-2022.03"
sbatch --mem=31G -n 24 --array=1-"$NUM_P19_Gut" $MC/scripts/gather_SLURM.sh "P19_Gut" $FECES_TSV "gtdb-rs214-rep"
sbatch --mem=80G -n 16 --array=1-"$NUM_P19_Gut" $MC/scripts/gather_SLURM.sh "P19_Gut" $FECES_TSV "gtdb-rs214-full"

sbatch --mem=80G -n 16 --array=1-"$NUM_Moss" $MC/scripts/gather_SLURM.sh "Moss" $MOSS_TSV "genbank-2022.03"
sbatch --mem=31G -n 24 --array=1-"$NUM_Moss" $MC/scripts/gather_SLURM.sh "Moss" $MOSS_TSV "gtdb-rs214-rep"
sbatch --mem=80G -n 16 --array=1-"$NUM_Moss" $MC/scripts/gather_SLURM.sh "Moss" $MOSS_TSV "gtdb-rs214-full"

sbatch --mem=120G -n 24 --array=1-"$NUM_NAFLD" $MC/scripts/gather_SLURM_fast.sh "NAFLD" $NAFLD_TSV "genbank-2022.03"
sbatch --mem=31G -n 24 --array=1-"$NUM_NAFLD" $MC/scripts/gather_SLURM_fast.sh "NAFLD" $NAFLD_TSV "gtdb-rs214-rep"
sbatch --mem=31G -n 24 --array=1-"$NUM_NAFLD" $MC/scripts/gather_SLURM_fast.sh "NAFLD" $NAFLD_TSV "gtdb-rs214-full"

sbatch --mem=120G -n 24--array=1-"$NUM_AD_Skin" $MC/scripts/gather_SLURM_fast.sh "AD_Skin" $AD_Skin_TSV "genbank-2022.03"
sbatch --mem=31G -n 24 --array=1-"$NUM_AD_Skin" $MC/scripts/gather_SLURM_fast.sh "AD_Skin" $AD_Skin_TSV "gtdb-rs214-rep"
sbatch --mem=80G -n 16 --array=1-"$NUM_AD_Skin" $MC/scripts/gather_SLURM_fast.sh "AD_Skin" $AD_Skin_TSV "gtdb-rs214-full"

sbatch --mem=120G -n 24 --array=1-"$NUM_PD" $MC/scripts/gather_SLURM_fast.sh "PD" $PD_TSV "genbank-2022.03"
sbatch --mem=31G -n 24 --array=1-"$NUM_PD" $MC/scripts/gather_SLURM_fast.sh "PD" $PD_TSV "gtdb-rs214-rep"
sbatch --mem=80G -n 16 --array=1-"$NUM_PD" $MC/scripts/gather_SLURM_fast.sh "PD" $PD_TSV "gtdb-rs214-full"

# Check completion status
check_output 'gtdb-rs214-rep gtdb-rs214-full genbank-2022.03' 'PD' _gather.csv

# Extract the lineage subset 
for i in $DATASETS; do
for j in SM_gtdb-rs214-full SM_gtdb-rs214-rep SM_genbank-2022.03; do
	db=$(echo "$j" | cut -d'_' -f2 | cut -d'-' -f1,2)
	ver=$(echo "$j" | cut -d'_' -f3)
cat $i/$j/*_gather.csv | cut -d, -f10 | tail -n+2 | awk '{print $1}' | sed 's/"//' | sort -u | \
	grep -Fhf - $ILAFORES/ref_dbs/sourmash_db/${db}*${ver}*.lineages.csv > $i/$j/${j}_lineages.csv
done
done

## surplus taxa in sourmash rs220 index
cat $MC/P19_Gut/Sourmash/*rs220*_gather.csv | cut -d, -f10 | tail -n+2 | \
	awk '{print $1}' | sed 's/"//' | sort -u > found_taxa.tsv # 8070 taxa

# of which 1191 are not in the species reps lineage file
grep -v -f <(cut -f1 $ILAFORES/ref_dbs/sourmash_db/bac120_taxonomy_r220.tsv | sed 's/^[^_]*_//' | sort -u) found_taxa.tsv | wc






## Remove and Redownload corrupted samples:
remove_these=($(sed -n "$(echo $missing_samples | sed 's/,/p;/g' | sed 's/;$//')" $dataset/preproc/preprocessed_reads.sample.tsv | awk '{print $1}' ))
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

