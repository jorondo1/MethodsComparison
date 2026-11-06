# Horse (PRJNA1102860) on Fir: raw reads to profiles. Run setup_fir.sh first.
# Copy-paste collection: running or sourcing this file does nothing. Paste a block, wait for it to finish, then the next.
return 0 2>/dev/null; [[ $0 == "${BASH_SOURCE[0]}" ]] && { echo "Copy-paste collection: paste the commands, don't run this file"; exit 1; }

source ~/scratch/MethodsComparison/scripts/fir/fir_env.sh
cd $MC && source scripts/src/myFunctions.sh

SRC=$MC/scripts/src
FIR=$MC/scripts/fir
DATASET_PATH=$MC/data/Horse
PREPROC=$DATASET_PATH/preproc
HOST_DB=$FAST/host_genomes/EquCab3_GRCh38_index/EquCab3_GRCh38
HORSE_INDEX=              # bowtie2 prefix of the horse index, under $FAST/host_genomes
DB_MPA=$FAST/metaphlan4_db
DB_SM=$ILAFORES/ref_dbs/sourmash_db
DB_SYLPH=$ILAFORES/ref_dbs/sylph_db
DB_K2=$ILAFORES/ref_dbs/kraken2_dbs

mkdir -p $DATASET_PATH/raw $PREPROC/logs $MC/logs

# run <job name> <sbatch args>: the tool scripts hardcode their -D and -o for the HPC, so they are set here
run() {
	local name=$1; shift
	sbatch -D $MC -J $name --mail-type=NONE -o $MC/logs/$name-%A_%a.out -e $MC/logs/$name-%A_%a.err "$@"
}

################################################################
# Sample sheet ########
################################################################

for r1 in $RAW/*_1.fastq.gz; do
	r2=${r1%_1.fastq.gz}_2.fastq.gz
	[[ -s $r2 ]] || { echo "Skipped, no pair: $r1" >&2; continue; }
	printf '%s\t%s\t%s\n' "$(basename "$r1" _1.fastq.gz)" "$r1" "$r2"
done > $DATASET_PATH/raw/samples_to_process.tsv
wc -l $DATASET_PATH/raw/samples_to_process.tsv

################################################################
# Host index: horse + human ########
################################################################

run host_index $FIR/build_host_index.sh $FAST/host_genomes/GRCh38_index/grch38_1kgmaj $FAST/host_genomes/$HORSE_INDEX $HOST_DB

################################################################
# Kneaddata ########
################################################################

bash $ILL_PIPELINES/generateslurm_preprocess.kneaddata.sh \
	--sample_tsv $DATASET_PATH/raw/samples_to_process.tsv \
	--out $PREPROC \
	--trimmomatic_options "SLIDINGWINDOW:4:20 MINLEN:50" \
	--db $HOST_DB \
	--slurm_mem 40G --slurm_threads 8

N_RAW=$(wc -l < $DATASET_PATH/raw/samples_to_process.tsv)
sbatch --time=6:00:00 --array=1-$N_RAW $PREPROC/preprocess.kneaddata.slurm.sh

# Once finished, check a kneaddata log and the output sizes, then write the sample tsv for the tools:
while read -r s _; do
	p=$PREPROC/$s/$s
	[[ -s ${p}_paired_1.fastq.gz ]] || { echo "Missing kneaddata output: $s" >&2; continue; }
	printf '%s\t%s\t%s\t%s\t%s\n' "$s" "${p}_paired_1.fastq.gz" "${p}_paired_2.fastq.gz" "${p}_unmatched_1.fastq.gz" "${p}_unmatched_2.fastq.gz"
done < $DATASET_PATH/raw/samples_to_process.tsv > $PREPROC/preprocessed_reads.sample.tsv

use_dataset Horse
CHUNK=1-$N_SAMPLES # if the job limit refuses an array, submit it in parts: CHUNK=1-400, then CHUNK=401-$N_SAMPLES

check_all Horse    # progress, between the blocks below
check_all -m Horse # + the unfinished rows, for CHUNK

################################################################
# One block at a time: each array is N_SAMPLES queued jobs ########
################################################################

# mOTUs ----------------
run motus4 --array=$CHUNK -n 8 --mem=80G --time=6:00:00 $SRC/motus4_SLURM.sh $DATASET_PATH $TSV
run motus3 --array=$CHUNK -n 8 --mem=40G --time=6:00:00 $SRC/motus3_SLURM.sh $DATASET_PATH $TSV

# MetaPhlAn ----------------
metaphlan="bash $ILL_PIPELINES/generateslurm_taxonomic_abundance.metaphlan.4.2.4.sh \
	--slurm_log $MC/logs --slurm_walltime 12:00:00 --slurm_threads 8 --slurm_mem 80G --sample_tsv $TSV"
$metaphlan --db $DB_MPA/2022/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $DATASET_PATH/MPA_db2022
sbatch --array=$CHUNK $DATASET_PATH/MPA_db2022/metaphlan.slurm.sh

$metaphlan --db $DB_MPA/2023/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $DATASET_PATH/MPA_db2023
sbatch --array=$CHUNK $DATASET_PATH/MPA_db2023/metaphlan.slurm.sh

$metaphlan --db $DB_MPA/2025/mpa_vJan25_CHOCOPhlAnSGB_202503 --out $DATASET_PATH/MPA_db2025
sbatch --array=$CHUNK $DATASET_PATH/MPA_db2025/metaphlan.slurm.sh

# sylph ----------------
run sylph --array=$CHUNK -n 4 --mem=64G --time=3:00:00 $SRC/sylph_SLURM.sh \
	Horse $TSV $DB_SYLPH/RefSeq_20250528_c200.syldb $DB_SYLPH/RefSeq_20250528_taxonomy.tsv.gz RefSeq
run sylph --array=$CHUNK -n 4 --mem=64G --time=3:00:00 $SRC/sylph_SLURM.sh \
	Horse $TSV $DB_SYLPH/gtdb-r220-c200-dbv1.syldb $DB_SYLPH/gtdb_r220_metadata.tsv.gz gtdb_rs220

# sourmash ----------------
# Pilot first: this script overwrote a read file once, on the first task of a new dataset. Fixed, but run CHUNK=1 alone,
# and check that it wrote data/Horse/signatures/ and left the reads alone.
run sourmash --array=$CHUNK -n 8 --mem=64G --time=12:00:00 $SRC/gather_SLURM_fast.sh Horse $TSV $DB_SM/gtdb-rs220-reps-k31.zip
run sourmash --array=$CHUNK -n 1 --mem=64G --time=12:00:00 $SRC/gather_SLURM_fast.sh Horse $TSV $DB_SM/RefSeq_20250528.k31.sbt.zip

# Kraken ----------------
run kraken $FIR/kraken_node.sh $DB_K2/k2_standard_20241228 "" $TSV $DATASET_PATH
run kraken $FIR/kraken_node.sh $DB_K2/k2_gtdb_genome_reps_20241109 _GTDB $TSV $DATASET_PATH

################################################################
# Results ########
################################################################

check_all Horse
run merge_tables $SRC/merge_tables.sh Horse

# To project, without reads and alignments
mkdir -p $RESULTS
rsync -a --no-g --no-p --exclude preproc --exclude raw --exclude signatures --exclude sylph_sketches_c200 \
	--exclude '*.bam' --exclude '*.inserts.gz' --exclude '*.tmp.csv' $DATASET_PATH/ $RESULTS/
