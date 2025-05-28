ml apptainer

sourmash="singularity exec --writable-tmpfs -e -B $ILAFORES:$ILAFORES,/fast2/def-ilafores:/fast2/def-ilafores $ILL_PIPELINES/containers/sourmash.4.8.11.sif sourmash"
refseq_genomes="/fast2/def-ilafores/refseq_genomes"
dwnld_date=refseq-09-April-2025

find $PWD/split_genomes_09-April-2025/signatures/ -type f -name "*.sig" > split_genomes_09-April-2025/signature_list.txt

singularity exec --writable-tmpfs -e -B $ANCHOR$ILAFORES:$ANCHOR$ILAFORES,$ANCHOR/fast2/def-ilafores:$ANCHOR/fast2/def-ilafores $ANCHOR$ILL_PIPELINES/containers/sourmash.4.8.11.sif sourmash index -k=31 split_genomes_09-April-2025/RefSeq_abpv_09042025 --from-file $ANCHOR$PWD/split_genomes_09-April-2025/signature_list.txt


# branchwater multithread sketch

# Build csv
echo name,genome_filename,protein_filename > $refseq_genomes/manysketch_${dwnld_date}.csv
for i in ${refseq_genomes}/*/${dwnld_date}-*-genomes/*.genomic.fna.gz
do
echo $i,$i,
done >> $refseq_genomes/manysketch_${dwnld_date}.csv

sourmash scripts manysketch fa.csv -o fa.zip -p k=21,k=31,k=51,scaled=1000,abund -p protein,k=10,scaled=200

$sourmash scripts manysketch $refseq_genomes/manysketch_${dwnld_date}.csv \
--cores 16 -p k=31,scaled=1000,abund --singleton -o $refseq_genomes/all_sig_${dwnld_date}.zip

######### OLD SINGLE THREAD
# Compute signature 
for bibitte in archaea bacteria plasmid viral; do
	mkdir -p ${refseq_genomes}/${bibitte}/${dwnld_date}-signatures
	$sourmash sketch dna ${refseq_genomes}/${bibitte}/${dwnld_date}-${bibitte}-genomes/*.genomic.fna.gz \
	-p k=31,scaled=1000,abund \
	--output-dir ${refseq_genomes}/${bibitte}/${dwnld_date}-signatures \
	--singleton # >> singleton means one signature for each record
done

# Create index
$sourmash sig cat $refseq_genomes/*/${dwnld_date}-signatures/*.sig -o $refseq_genomes/all_sig.zip
# OR index (may take weeks)
$sourmash index -k 31 $refseq_genomes/RefSeq-ABPV-090425.k31.sbt.zip $refseq_genomes/*/${dwnld_date}-signatures/*.sig

cp $refseq_genomes/RefSeq-ABPV-090425.k31.rocksdb

# Extract taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.EXTRA.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.Z

uncompress -c taxdump.tar.Z | tar xf - 

