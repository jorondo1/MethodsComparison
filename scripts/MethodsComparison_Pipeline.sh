export MC=$ILAFORES/analysis/MethodsComparison
export ILL_PIPELINES=$ILAFORES/analysis/MethodsComparison/ILL_pipelines
cd $MC

export SALIVA_TSV=$PR19/Saliva/preproc/preprocessed_reads.sample.tsv
export FECES_TSV=$PR19/Feces/preproc/preprocessed_reads.sample.tsv
export MOSS_TSV=$MOSS/preproc/preprocessed_reads.sample.tsv
export NUM_Saliva=$(wc $SALIVA_TSV | awk '{print $1}')
export NUM_Feces=$(wc $FECES_TSV | awk '{print $1}')
export NUM_Moss=$(wc $MOSS_TSV | awk '{print $1}')
export DATASETS="Saliva Feces Moss"

############
# mOTUs ####
# Custom SLURM script 
sbatch --array=1-"$NUM_Saliva" $MC/scripts/motus_SLURM.sh Saliva $SALIVA_TSV
sbatch --array=1-"$NUM_Feces" $MC/scripts/motus_SLURM.sh Feces $FECES_TSV
sbatch --array=1-"$NUM_Moss" $MC/scripts/motus_SLURM.sh Moss $MOSS_TSV

# Check completion status
for i in $DATASETS; do
	eval exp=\$NUM_$i
	num=$(ls $i/mOTU_abund/*_profile.txt | wc -l)
	echo "$num mOTUs output for $i found, $exp expected."
done

####################
# Kraken/bracken ###
bracken="bash $ILL_PIPELINES/generateslurm_taxonomic_profile.sample.sh \
	--kraken_db $ILAFORES/ref_dbs/kraken2_dbs/kraken2_PlusPFP_202202 \
	--slurm_log $MC/logs --slurm_walltime 72:00:00 --slurm_threads 48 --slurm_mem 250G"

# Generate SLURM scripts https://github.com/jflucier/ILL_pipelines/blob/main/generateslurm_taxonomic_profile.sample.sh
$bracken --confidence 0.20 --sample_tsv $SALIVA_TSV --out $MC/Saliva/KB20
$bracken --confidence 0.20 --sample_tsv $FECES_TSV --out $MC/Feces/KB20
$bracken --confidence 0.20 --sample_tsv $MOSS_TSV --out $MC/Moss/KB20

$bracken --confidence 0.51 --sample_tsv $SALIVA_TSV --out $MC/Saliva/KB51
$bracken --confidence 0.51 --sample_tsv $FECES_TSV --out $MC/Feces/KB51
$bracken --confidence 0.51 --sample_tsv $MOSS_TSV --out $MC/Moss/KB51

# Submit
sbatch --array=1-"$NUM_Saliva" $MC/Saliva/KB20/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_Feces" $MC/Feces/KB20/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_Moss" $MC/Moss/KB20/taxonomic_profile.samples.slurm.sh

sbatch --array=1-"$NUM_Saliva" $MC/Saliva/KB51/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_Feces" $MC/Feces/KB51/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_Moss" $MC/Moss/KB51/taxonomic_profile.samples.slurm.sh

# Check completion status
for test in KB20; do
for i in $DATASETS; do
	eval exp=\$NUM_$i
	num=$(ls $i/$test/*/*_bracken/*_bracken_S.MPA.TXT | wc -l)
	echo "$num $test output for $i found, $exp expected."
done
done

# Once completely done, remove taxonomy_nt files from kraken out 
rm */Bracken*/*/*_taxonomy_nt 
rm */Bracken*/*/*/*.bracken
rm */Bracken*/*/*/*.kreport
rm */Bracken*/*/*bugs_list.MPA.TXT
rm */Bracken*/*/*_temp.MPA.TXT

################
# MetaPhlAn4 ###
metaphlan="bash $ILL_PIPELINES/generateslurm_taxonomic_abundance.metaphlan.sh \
	--slurm_log $MC/logs --slurm_walltime 24:00:00 --slurm_threads 24 --slurm_mem 30G"

# Generate SLURM scripts https://github.com/jflucier/ILL_pipelines/blob/main/generateslurm_taxonomic_abundance.metaphlan.sh
# 2022 database
$metaphlan --sample_tsv $SALIVA_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/Saliva/MPA_db2022
$metaphlan --sample_tsv $FECES_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/Feces/MPA_db2022
$metaphlan --sample_tsv $MOSS_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/Moss/MPA_db2022

sbatch --array=1-"$NUM_Saliva" $MC/Saliva/MPA_db2022/metaphlan.slurm.sh
sbatch --array=1-"$NUM_Feces" $MC/Feces/MPA_db2022/metaphlan.slurm.sh
sbatch --array=1-"$NUM_Moss" $MC/Moss/MPA_db2022/metaphlan.slurm.sh


# 2023 database
$metaphlan --sample_tsv $SALIVA_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/Saliva/MPA_db2023
$metaphlan --sample_tsv $FECES_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/Feces/MPA_db2023
$metaphlan --sample_tsv $MOSS_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/Moss/MPA_db2023

sbatch --array=1-"$NUM_Saliva" $MC/Saliva/MPA_db2023/metaphlan.slurm.sh
sbatch --array=1-"$NUM_Feces" $MC/Feces/MPA_db2023/metaphlan.slurm.sh
sbatch --array=1-"$NUM_Moss" $MC/Moss/MPA_db2023/metaphlan.slurm.sh

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

#####################
# Sourmash gather ###

sbatch --mem=60G --array=1-"$NUM_Saliva" $MC/scripts/gather_SLURM.sh "Saliva" $SALIVA_TSV "genbank-2022.03"
sbatch --mem=80G --array=1-"$NUM_Feces" $MC/scripts/gather_SLURM.sh "Feces" $FECES_TSV "genbank-2022.03"
sbatch --mem=80G --array=1-"$NUM_Moss" $MC/scripts/gather_SLURM.sh "Moss" $MOSS_TSV "genbank-2022.03"

sbatch --mem=31G --array=1-"$NUM_Saliva" $MC/scripts/gather_SLURM.sh "Saliva" $SALIVA_TSV "gtdb-rs214-rep"
sbatch --mem=31G --array=1-"$NUM_Feces" $MC/scripts/gather_SLURM.sh "Feces" $FECES_TSV "gtdb-rs214-rep"
sbatch --mem=31G --array=1-"$NUM_Moss" $MC/scripts/gather_SLURM.sh "Moss" $MOSS_TSV "gtdb-rs214-rep"

sbatch --mem=80G --array=1-"$NUM_Saliva" $MC/scripts/gather_SLURM.sh "Saliva" $SALIVA_TSV "gtdb-rs214-full"
sbatch --mem=80G --array=1-"$NUM_Feces" $MC/scripts/gather_SLURM.sh "Feces" $FECES_TSV "gtdb-rs214-full"
sbatch --mem=80G --array=1-"$NUM_Moss" $MC/scripts/gather_SLURM.sh "Moss" $MOSS_TSV "gtdb-rs214-full"

# Check completion status
for SM_db in gtdb_rs214_rep gtdb_rs214_full; do
for i in $DATASETS; do
	eval exp=\$NUM_$i
	num=$(ls $i/SM_${SM_db}/*${SM_db}_gather.csv | wc -l)
	echo "$num $SM_db output for $i found, $exp expected."
done
done

# Extract the lineage subset 
for i in $DATASETS; do
for j in SM_gtdb_rs214_full SM_gtdb_rs214_rep; do
	db=$(echo "$j" | cut -d'_' -f2)
	ver=$(echo "$j" | cut -d'_' -f3)
cat $i/$j/*_gather.csv | cut -d, -f10 | tail -n+2 | awk '{print $1}' | sed 's/"//' | sort -u | \
	grep -Fhf - $ILAFORES/ref_dbs/sourmash_db/${db}*${ver}*.lineages.csv > $i/$j/${j}_lineages.csv
done
done


## surplus taxa in sourmash rs220 index
cat $MC/Feces/Sourmash/*rs220*_gather.csv | cut -d, -f10 | tail -n+2 | \
	awk '{print $1}' | sed 's/"//' | sort -u > found_taxa.tsv # 8070 taxa

# of which 1191 are not in the species reps lineage file
grep -v -f <(cut -f1 $ILAFORES/ref_dbs/sourmash_db/bac120_taxonomy_r220.tsv | sed 's/^[^_]*_//' | sort -u) found_taxa.tsv | wc





# Process fecal samples 
# bash $ILL_PIPELINES/generateslurm_preprocess.kneaddata.sh \
# 	--sample_tsv $MC/Feces/samples_to_process.tsv \
# 	--out $ILAFORES/analysis/projet_PROVID19/Feces/preproc \
# 	--db $FAST/host_genomes/GRCh38_index/grch38_1kgmaj \
# 	--slurm_mem 125G --slurm_threads 24
# correct script as the shell command needs to be anchored!
