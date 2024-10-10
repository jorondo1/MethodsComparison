export MC=$ILAFORES/analysis/MethodsComparison
export ILL_PIPELINES=$ILAFORES/analysis/MethodsComparison/ILL_pipelines
cd $MC

export SALIVA_TSV=$PR19/P19_Saliva/preproc/preprocessed_reads.sample.tsv
export FECES_TSV=$PR19/P19_Gut/preproc/preprocessed_reads.sample.tsv
export MOSS_TSV=$MOSS/preproc/preprocessed_reads.sample.tsv
export NAFLD_TSV=$MC/NAFLD/preproc/preprocessed_reads.sample.tsv
export AD_Skin_TSV=$MC/AD_Skin/preproc/preprocessed_reads.sample.tsv
export NUM_P19_Saliva=$(wc $SALIVA_TSV | awk '{print $1}')
export NUM_P19_Gut=$(wc $FECES_TSV | awk '{print $1}')
export NUM_Moss=$(wc $MOSS_TSV | awk '{print $1}')
export NUM_NAFLD=$(wc $NAFLD_TSV | awk '{print $1}')
export NUM_AD_Skin=$(wc $AD_Skin_TSV | awk '{print $1}')
export DATASETS="P19_Saliva P19_Gut Moss NAFLD AD_Skin"


############
# mOTUs ####
# Custom SLURM script 
sbatch --array=1-"$NUM_P19_Saliva" $MC/scripts/motus_SLURM.sh P19_Saliva $SALIVA_TSV
sbatch --array=1-"$NUM_P19_Gut" $MC/scripts/motus_SLURM.sh P19_Gut $FECES_TSV
sbatch --array=1-"$NUM_Moss" $MC/scripts/motus_SLURM.sh Moss $MOSS_TSV
sbatch --array=1-"$NUM_NAFLD" $MC/scripts/motus_SLURM.sh NAFLD $NAFLD_TSV
sbatch --array=1-"$NUM_AD_Skin" $MC/scripts/motus_SLURM.sh NAFLD $AD_Skin_TSV

# Check completion status
for i in $DATASETS; do
	eval exp=\$NUM_$i # find which files have more than 1 line : 
	num=0
	for file in "$i/MOTUS"/*_profile.txt; do
	  if [ $(wc -l < "$file") -gt 1 ]; then
	    num=$((num + 1))
	  fi
	done
	echo "$num mOTUs output for $i found, $exp expected."
done

####################
# Kraken/bracken ###
bracken="bash $ILL_PIPELINES/generateslurm_taxonomic_profile.sample.sh \
	--kraken_db $ILAFORES/ref_dbs/kraken2_dbs/kraken2_PlusPFP_202202 \
	--slurm_log $MC/logs --slurm_walltime 72:00:00 --slurm_threads 48 --slurm_mem 250G"

# Generate SLURM scripts https://github.com/jflucier/ILL_pipelines/blob/main/generateslurm_taxonomic_profile.sample.sh
$bracken --confidence 0.05 --sample_tsv $SALIVA_TSV --out $MC/P19_Saliva/KB20
$bracken --confidence 0.05 --sample_tsv $FECES_TSV --out $MC/P19_Gut/KB20
$bracken --confidence 0.05 --sample_tsv $MOSS_TSV --out $MC/Moss/KB20
$bracken --confidence 0.05 --sample_tsv $NAFLD_TSV --out $MC/NAFLD/KB05

$bracken --confidence 0.20 --sample_tsv $SALIVA_TSV --out $MC/P19_Saliva/KB20
$bracken --confidence 0.20 --sample_tsv $FECES_TSV --out $MC/P19_Gut/KB20
$bracken --confidence 0.20 --sample_tsv $MOSS_TSV --out $MC/Moss/KB20
$bracken --confidence 0.20 --sample_tsv $NAFLD_TSV --out $MC/NAFLD/KB20
$bracken --confidence 0.20 --sample_tsv $AD_Skin_TSV --out $MC/AD_Skin/KB20

$bracken --confidence 0.51 --sample_tsv $SALIVA_TSV --out $MC/P19_Saliva/KB51
$bracken --confidence 0.51 --sample_tsv $FECES_TSV --out $MC/P19_Gut/KB51
$bracken --confidence 0.51 --sample_tsv $MOSS_TSV --out $MC/Moss/KB51
$bracken --confidence 0.51 --sample_tsv $NAFLD_TSV --out $MC/NAFLD/KB51
$bracken --confidence 0.51 --sample_tsv $AD_Skin_TSV --out $MC/AD_Skin/KB51

# Submit
sbatch --array=1-"$NUM_P19_Saliva" $MC/P19_Saliva/KB05/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_P19_Gut" $MC/P19_Gut/KB05/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_Moss" $MC/Moss/KB05/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_NAFLD" $MC/NAFLD/KB05/taxonomic_profile.samples.slurm.sh

sbatch --array=1-"$NUM_P19_Saliva" $MC/P19_Saliva/KB20/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_P19_Gut" $MC/P19_Gut/KB20/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_Moss" $MC/Moss/KB20/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_NAFLD" $MC/NAFLD/KB20/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_AD_Skin" $MC/AD_Skin/KB20/taxonomic_profile.samples.slurm.sh

sbatch --array=1-"$NUM_P19_Saliva" $MC/P19_Saliva/KB51/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_P19_Gut" $MC/P19_Gut/KB51/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_Moss" $MC/Moss/KB51/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_NAFLD" $MC/NAFLD/KB51/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_AD_Skin" $MC/AD_Skin/KB51/taxonomic_profile.samples.slurm.sh

# Check completion status
for test in KB20 KB51; do
for i in $DATASETS; do
	eval exp=\$NUM_$i
	num=$(ls $i/$test/*/*_bracken/*_bracken_S.MPA.TXT | wc -l)
	echo "$num $test output for $i found, $exp expected."
done
done

# Once completely done, remove taxonomy_nt files from kraken out 
rm */KB*/*/*_taxonomy_nt 
rm */KB*/*/*/*.bracken
rm */KB*/*/*/*.kreport
rm */KB*/*/*bugs_list.MPA.TXT
rm */KB*/*/*/*_temp.MPA.TXT

################
# MetaPhlAn4 ###
metaphlan="bash $ILL_PIPELINES/generateslurm_taxonomic_abundance.metaphlan.sh \
	--slurm_log $MC/logs --slurm_walltime 24:00:00 --slurm_threads 24 --slurm_mem 30G"

# Generate SLURM scripts https://github.com/jflucier/ILL_pipelines/blob/main/generateslurm_taxonomic_abundance.metaphlan.sh
# 2022 database
$metaphlan --sample_tsv $SALIVA_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/P19_Saliva/MPA_db2022
$metaphlan --sample_tsv $FECES_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/P19_Gut/MPA_db2022
$metaphlan --sample_tsv $MOSS_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/Moss/MPA_db2022
$metaphlan --sample_tsv $NAFLD_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/NAFLD/MPA_db2022
$metaphlan --sample_tsv $AD_Skin_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/AD_Skin/MPA_db2022

sbatch --array=1-"$NUM_P19_Saliva" $MC/P19_Saliva/MPA_db2022/metaphlan.slurm.sh
sbatch --array=1-"$NUM_P19_Gut" $MC/P19_Gut/MPA_db2022/metaphlan.slurm.sh
sbatch --array=1-"$NUM_Moss" $MC/Moss/MPA_db2022/mfetaphlan.slurm.sh
sbatch --array=1-"$NUM_NAFLD" $MC/NAFLD/MPA_db2022/metaphlan.slurm.sh
sbatch --array=1-"$NUM_AD_Skin" $MC/AD_Skin/MPA_db2022/metaphlan.slurm.sh

# 2023 database
$metaphlan --sample_tsv $SALIVA_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/P19_Saliva/MPA_db2023
$metaphlan --sample_tsv $FECES_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/P19_Gut/MPA_db2023
$metaphlan --sample_tsv $MOSS_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/Moss/MPA_db2023
$metaphlan --sample_tsv $NAFLD_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/NAFLD/MPA_db2023
$metaphlan --sample_tsv $AD_Skin_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/AD_Skin/MPA_db2023

sbatch --array=1-"$NUM_P19_Saliva" $MC/P19_Saliva/MPA_db2023/metaphlan.slurm.sh
sbatch --array=1-"$NUM_P19_Gut" $MC/P19_Gut/MPA_db2023/metaphlan.slurm.sh
sbatch --array=1-"$NUM_Moss" $MC/Moss/MPA_db2023/metaphlan.slurm.sh
sbatch --array=1-"$NUM_NAFLD" $MC/NAFLD/MPA_db2023/metaphlan.slurm.sh
sbatch --array=1-"$NUM_AD_Skin" $MC/AD_Skin/MPA_db2023/metaphlan.slurm.sh

# Check completion status
for MPA in MPA_db2022 MPA_db2023; do
for i in $DATASETS; do
	eval exp=\$NUM_$i
	num=$(ls $i/$MPA/*/*_profile.txt | wc -l)
	echo "$num $MPA output for $i found, $exp expected."
done
done

# Remove bowtie indexes
rm */MPA_db*/*/*.bowtie2.txt
rm */*/.throttle -r

#####################
# Sourmash gather ###

sbatch --mem=60G --array=1-"$NUM_P19_Saliva" $MC/scripts/gather_SLURM.sh "P19_Saliva" $SALIVA_TSV "genbank-2022.03"
sbatch --mem=80G --array=1-"$NUM_P19_Gut" $MC/scripts/gather_SLURM.sh "P19_Gut" $FECES_TSV "genbank-2022.03"
sbatch --mem=80G --array=1-"$NUM_Moss" $MC/scripts/gather_SLURM.sh "Moss" $MOSS_TSV "genbank-2022.03"
sbatch --mem=120G --array=1-"$NUM_NAFLD" $MC/scripts/gather_SLURM.sh "NAFLD" $NAFLD_TSV "genbank-2022.03"
sbatch --mem=120G --array=1-"$NUM_AD_Skin" $MC/scripts/gather_SLURM.sh "AD_Skin" $NAFLD_TSV "genbank-2022.03"

sbatch --mem=31G --array=1-"$NUM_P19_Saliva" $MC/scripts/gather_SLURM.sh "P19_Saliva" $SALIVA_TSV "gtdb-rs214-rep"
sbatch --mem=31G --array=1-"$NUM_P19_Gut" $MC/scripts/gather_SLURM.sh "P19_Gut" $FECES_TSV "gtdb-rs214-rep"
sbatch --mem=31G --array=1-"$NUM_Moss" $MC/scripts/gather_SLURM.sh "Moss" $MOSS_TSV "gtdb-rs214-rep"
sbatch --mem=31G --array=1-"$NUM_NAFLD" $MC/scripts/gather_SLURM.sh "NAFLD" $NAFLD_TSV "gtdb-rs214-rep"
sbatch --mem=31G --array=1-"$NUM_AD_Skin" $MC/scripts/gather_SLURM.sh "AD_Skin" $AD_Skin_TSV "gtdb-rs214-rep"

sbatch --mem=80G --array=1-"$NUM_P19_Saliva" $MC/scripts/gather_SLURM.sh "P19_Saliva" $SALIVA_TSV "gtdb-rs214-full"
sbatch --mem=80G --array=1-"$NUM_P19_Gut" $MC/scripts/gather_SLURM.sh "P19_Gut" $FECES_TSV "gtdb-rs214-full"
sbatch --mem=80G --array=1-"$NUM_Moss" $MC/scripts/gather_SLURM.sh "Moss" $MOSS_TSV "gtdb-rs214-full"
sbatch --mem=31G --array=1-"$NUM_NAFLD" $MC/scripts/gather_SLURM.sh "NAFLD" $NAFLD_TSV "gtdb-rs214-full"
sbatch --mem=31G --array=1-"$NUM_AD_Skin" $MC/scripts/gather_SLURM.sh "AD_Skin" $AD_Skin_TSV "gtdb-rs214-full"

# Check completion status
for SM_db in gtdb_rs214_rep gtdb_rs214_full genbank-2022.03; do
for i in $DATASETS; do
	eval exp=\$NUM_$i
	num=$(ls $i/SM_*/*${SM_db}_gather.csv | wc -l)
	echo "$num $SM_db output for $i found, $exp expected."
done
done

# Extract the lineage subset 
for i in $DATASETS; do
for j in SM_gtdb_rs214_full SM_gtdb_rs214_rep SM_genbank-2022.03; do
	db=$(echo "$j" | cut -d'_' -f2)
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

######################
# generate preprocess file 

#Process fecal samples
mkdir -p $MC/NAFLD/preproc
bash $ILAFORES/programs/ILL_pipelines/generateslurm_preprocess.kneaddata.sh \
	--sample_tsv $MC/NAFLD/raw/samples_to_process.tsv \
	--out $MC/NAFLD/preproc \
	--trimmomatic_options "SLIDINGWINDOW:4:20 MINLEN:50" \
	--db $FAST/host_genomes/GRCh38_index/grch38_1kgmaj \
	--slurm_mem 31G --slurm_threads 24
# correct script as the shell command needs to be anchored!

# Find missing line numbers
grep -n -v -f <(ls NAFLD/SM_genbank_202203/ | sed 's/_.*//' | sort | uniq) NAFLD/preproc/preprocessed_reads.sample.tsv | awk -F: '{print $1}' | paste -sd,

:>$MC/NAFLD/raw/samples_to_process.tsv
while read -r p; do 
SRR=$(echo $p | awk '{print $1}')
f1=$(find NAFLD/raw -type f -name "${SRR}_1.fastq*" -exec realpath {} \;)
f2=$(find NAFLD/raw -name "${SRR}_2.fastq*" -exec realpath {} \;)
SAM=$(echo $p | awk '{print $13}')
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
