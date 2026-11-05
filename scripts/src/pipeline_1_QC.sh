# QC: download, kneaddata, gzip. Copy-paste collection: running or sourcing this file does nothing.
return 0 2>/dev/null; [[ $0 == "${BASH_SOURCE[0]}" ]] && { echo "Copy-paste collection: paste the commands, don't run this file"; exit 1; }

export MC=$ILAFORES/analysis/MethodsComparison
export MAPLE_DATA=/net/nfs-bio/jbod2/def-ilafores/analysis/mapleRhizosphere/WMS
cd $MC && source scripts/src/myFunctions.sh
use_dataset Maple # names: scripts/src/datasets.tsv (add the new dataset there first)

PARTS=bio,c-3t # + cpubase for light jobs (old nodes)
PREPROC=$DATASET_PATH/preproc # kneaddata output, where datasets.tsv points

################################################################################################################
# Public dataset (ENA) #################
################################################################################################################

mkdir -p $DATASET_PATH/raw $FAST/$DATASET/raw

# Download from ENA tsv file with grep subset ------
# ENA_DL_SCRIPT: the download script from the ENA browser ("Download all", bash script), saved in raw/
samples=($(grep "WGS" $DATASET_PATH/raw/filereport_*_tsv.txt | cut -f1 | grep -v 'run'))
ENA_DL_SCRIPT=$(ls $DATASET_PATH/raw/ena-file-download-*.sh)
cd $FAST/$DATASET/raw
grep -E "$(IFS="|"; echo "${samples[*]}")" $ENA_DL_SCRIPT | bash
cd $MC

# Create TSV
:> $DATASET_PATH/raw/samples_to_process.tsv
for sample in "${samples[@]}"; do
	fq1=$(find $FAST/$DATASET/raw -type f -name "${sample}_1.fastq.gz")
	fq2=$(find $FAST/$DATASET/raw -type f -name "${sample}_2.fastq.gz")
	echo -e "${sample}\t${fq1}\t${fq2}" >> $DATASET_PATH/raw/samples_to_process.tsv
done

# Process samples
mkdir -p $PREPROC/logs
bash $ILAFORES/programs/ILL_pipelines/generateslurm_preprocess.kneaddata.sh \
	--sample_tsv $DATASET_PATH/raw/samples_to_process.tsv \
	--out $PREPROC \
	--trimmomatic_options "SLIDINGWINDOW:4:20 MINLEN:50" \
	--db $FAST/host_genomes/GRCh38_index/grch38_1kgmaj \
	--slurm_mem 30G --slurm_threads 24

sbatch --array=1-$(wc -l < $DATASET_PATH/raw/samples_to_process.tsv) --partition $PARTS $PREPROC/preprocess.kneaddata.slurm.sh

# Find arrays of missing samples:
missing_samples=$(grep -n -v -wFf <(ls $PREPROC/*/*_1.fastq* | xargs -n1 dirname | xargs -n1 basename) $DATASET_PATH/raw/samples_to_process.tsv | cut -f1 -d: | paste -sd,); echo $missing_samples

rm -r $PREPROC/.throttle
sbatch --array="$missing_samples" --partition $PARTS $PREPROC/preprocess.kneaddata.slurm.sh

# gzip all, except if gz already there
find $PREPROC -type f -name '*fastq' -print0 | xargs -0 -P 8 -I {} bash -c 'if [ ! -f "{}.gz" ]; then pigz -k "{}"; fi'

use_dataset $DATASET # now that the sample tsv exists

# ── ONE-TIME (done): .fastq -> .fastq.gz and path fixes in every sample tsv ──
for tsv in $(find $ILAFORES/analysis/ -maxdepth 5 -readable -name "preprocessed_reads.sample.tsv" 2>/dev/null); do
sed -i 's|.fastq.gz|.fastq|g' $tsv
sed -i 's|.fastq|.fastq.gz|g' $tsv
sed -i 's|/nfs3_ib/nfs-ip34/home|/net/nfs-ip34/jbod2|g' $tsv
sed -i 's|/net/nfs-ip34/home|/net/nfs-ip34/jbod2|g' $tsv
done

################################################################################################################
# Local not preprocessed dataset (Maple) #################
################################################################################################################

mapfile -t samples < <(tail -n +2 ${MAPLE_DATA}/raw/OMX753-WGS-13-05-2025_IndexAssignment.csv | grep WGS - | cut -d, -f2 | sort -u)

:> $MAPLE_DATA/raw/samples_to_process.tsv
for sample in "${samples[@]}"; do
	fq1=$(find $MAPLE -type f -name "${sample}_R1.fastq.gz")
	fq2=$(find $MAPLE -type f -name "${sample}_R2.fastq.gz")
	echo -e "${sample}\t${fq1}\t${fq2}" >> $MAPLE_DATA/raw/samples_to_process.tsv
done

# ── ONE-TIME (done): Acer saccharum host genome ──
mkdir -p $FAST/host_genomes/_genome_fasta/AcerSaccharum
for chromosome in CM061734.1 CM061735.1 CM061736.1 CM061737.1 CM061738.1 CM061739.1 CM061740.1 CM061741.1 CM061742.1 CM061743.1 CM061744.1 CM061745.1 CM061746.1; do
	efetch -db nucleotide -id ${chromosome} -format fasta > $FAST/host_genomes/_genome_fasta/AcerSaccharum/${chromosome}.fasta
done
# index
ml bowtie2
cat $FAST/host_genomes/_genome_fasta/AcerSaccharum/*.fasta >  $FAST/host_genomes/_genome_fasta/AcerSaccharum/full_AcerSaccharum.fasta
bowtie2-build --threads 48 --seed 1 $FAST/host_genomes/_genome_fasta/AcerSaccharum/full_AcerSaccharum.fasta $FAST/host_genomes/AcerSaccharum_index/AcerSaccharum
# ────

mkdir -p $MAPLE_DATA/preproc
bash $ILAFORES/programs/ILL_pipelines/generateslurm_preprocess.kneaddata.sh \
	--sample_tsv $MAPLE_DATA/raw/samples_to_process.tsv \
	--out $MAPLE_DATA/preproc \
	--trimmomatic_options "SLIDINGWINDOW:4:20 MINLEN:50" \
	--db $FAST/host_genomes/AcerSaccharum_index/AcerSaccharum \
	--slurm_mem 60G --slurm_threads 24

sbatch --partition $PARTS,cpubase --array=1-88 $MAPLE_DATA/preproc/preprocess.kneaddata.slurm.sh

missing_samples=$(grep -n -v -wFf <(ls $MAPLE_DATA/preproc/*/*R1.fastq* | xargs -n1 dirname | xargs -n1 basename) $MAPLE_DATA/raw/samples_to_process.tsv | cut -f1 -d: | paste -sd,); echo $missing_samples
sbatch --partition $PARTS,cpubase --array="$missing_samples" $MAPLE_DATA/preproc/preprocess.kneaddata.slurm.sh

# gzip output
sed -i 's|.fastq|.fastq.gz|g' $MAPLE_DATA/preproc/preprocessed_reads.sample.tsv
# ── DESTRUCTIVE: removes the existing .gz before recompressing ──
rm $MAPLE_DATA/preproc/*/*.gz
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

use_dataset Maple
