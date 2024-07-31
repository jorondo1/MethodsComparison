export MC=$ILAFORES/analysis/MethodsComparison
export ILL_PIPELINES=$ILAFORES/analysis/MethodsComparison/ILL_pipelines
cd $MC

# Process fecal samples 
# bash $ILL_PIPELINES/generateslurm_preprocess.kneaddata.sh \
# 	--sample_tsv $MC/Feces/samples_to_process.tsv \
# 	--out $ILAFORES/analysis/projet_PROVID19/Feces/preproc \
# 	--db $FAST/host_genomes/GRCh38_index/grch38_1kgmaj \
# 	--slurm_mem 125G --slurm_threads 24
# correct script as the shell command needs to be anchored!

export SALIVA_TSV=$PR19/Saliva/preproc/preprocessed_reads.sample.tsv
export FECES_TSV=$PR19/Feces/preproc/preprocessed_reads.sample.tsv
export MOSS_TSV=$MOSS/preproc/preprocessed_reads.sample.tsv
export NUM_SALIVA=$(wc $SALIVA_TSV | awk '{print $1}')
export NUM_FECES=$(wc $FECES_TSV | awk '{print $1}')
export NUM_MOSS=$(wc $MOSS_TSV | awk '{print $1}')

############
# mOTUs ####
# Custom SLURM script 
sbatch --array=1-"$NUM_SALIVA" $MC/scripts/motus_SLURM.sh Saliva $SALIVA_TSV
sbatch --array=1-"$NUM_FECES" $MC/scripts/motus_SLURM.sh Feces $FECES_TSV
sbatch --array=1-"$NUM_MOSS" $MC/scripts/motus_SLURM.sh Moss $MOSS_TSV

ls */mOTU_abund/*_profile.txt | wc -l

####################
# Kraken/bracken ###
bracken="bash $ILL_PIPELINES/generateslurm_taxonomic_profile.sample.sh \
	--kraken_db $ILAFORES/ref_dbs/kraken2_dbs/kraken2_PlusPFP_202202 \
	--slurm_log $MC/logs --slurm_walltime 24:00:00 --slurm_threads 48 --slurm_mem 250G"

# Generate SLURM scripts https://github.com/jflucier/ILL_pipelines/blob/main/generateslurm_taxonomic_profile.sample.sh
$bracken --sample_tsv $SALIVA_TSV --out $MC/Saliva/Bracken
$bracken --sample_tsv $FECES_TSV --out $MC/Feces/Bracken
$bracken --sample_tsv $MOSS_TSV --out $MC/Moss/Bracken

# Submit
sbatch --array=1-"$NUM_SALIVA" $MC/Saliva/Bracken/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_FECES" $MC/Feces/Bracken/taxonomic_profile.samples.slurm.sh
sbatch --array=1-"$NUM_MOSS" $MC/Moss/Bracken/taxonomic_profile.samples.slurm.sh

ls */Bracken/*/*-bugs_list.MPA.TXT | wc -l

################
# MetaPhlAn4 ###
metaphlan="bash $ILL_PIPELINES/generateslurm_taxonomic_abundance.metaphlan.sh \
	--slurm_log $MC/logs --slurm_walltime 24:00:00 --slurm_threads 24 --slurm_mem 30G"

# Generate SLURM scripts https://github.com/jflucier/ILL_pipelines/blob/main/generateslurm_taxonomic_abundance.metaphlan.sh
# 2022 database
$metaphlan --sample_tsv $SALIVA_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/Saliva/MPA_db2022
$metaphlan --sample_tsv $FECES_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/Feces/MPA_db2022
$metaphlan --sample_tsv $MOSS_TSV --db $FAST/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $MC/Moss/MPA_db2022

sbatch --array=1-"$NUM_SALIVA" $MC/Saliva/MPA_db2022/metaphlan.slurm.sh
sbatch --array=1-"$NUM_FECES" $MC/Feces/MPA_db2022/metaphlan.slurm.sh
sbatch --array=1-"$NUM_MOSS" $MC/Moss/MPA_db2022/metaphlan.slurm.sh

ls */MPA_db2022/*/*_profile.txt | wc -l

# 2023 database
$metaphlan --sample_tsv $SALIVA_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/Saliva/MPA_db2023
$metaphlan --sample_tsv $FECES_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/Feces/MPA_db2023
$metaphlan --sample_tsv $MOSS_TSV --db $FAST/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $MC/Moss/MPA_db2023

sbatch --array=1-"$NUM_SALIVA" $MC/Saliva/MPA_db2023/metaphlan.slurm.sh
sbatch --array=1-"$NUM_FECES" $MC/Feces/MPA_db2023/metaphlan.slurm.sh
sbatch --array=1-"$NUM_MOSS" $MC/Moss/MPA_db2023/metaphlan.slurm.sh

ls */MPA_db2023/*/*_profile.txt | wc -l

#####################
# Sourmash gather ###

sbatch --mem=60G --array=1-"$NUM_SALIVA" $MC/scripts/gather_SLURM.sh Saliva $SALIVA_TSV "genbank-2022.03"
sbatch --mem=60G --array=1-"$NUM_FECES" $MC/scripts/gather_SLURM.sh Feces $FECES_TSV "genbank-2022.03"
sbatch --mem=80G --array=1-"$NUM_MOSS" $MC/scripts/gather_SLURM.sh Moss $MOSS_TSV "genbank-2022.03"

sbatch --mem=15G --array=1-"$NUM_SALIVA" $MC/scripts/gather_SLURM.sh Saliva $SALIVA_TSV "gtdb-rs220"
sbatch --mem=15G --array=1-"$NUM_FECES" $MC/scripts/gather_SLURM.sh Feces $FECES_TSV "gtdb-rs220"
sbatch --mem=15G --array=1-"$NUM_MOSS" $MC/scripts/gather_SLURM.sh Moss $MOSS_TSV "gtdb-rs220"

















