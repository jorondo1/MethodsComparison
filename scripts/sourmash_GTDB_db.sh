ml apptainer

sourmash="singularity exec --writable-tmpfs -e -B $ILAFORES:$ILAFORES,/fast2/def-ilafores:/fast2/def-ilafores $ILL_PIPELINES/containers/sourmash.4.8.11.sif sourmash"
refseq_genomes="/fast2/def-ilafores/refseq_genomes"
dwnld_date=refseq-09-April-2025

# Compute signature 
for bibitte in archaea bacteria plasmid viral; do
	mkdir -p ${refseq_genomes}/${bibitte}/${dwnld_date}-signatures
	$sourmash sketch dna ${refseq_genomes}/${bibitte}/${dwnld_date}-${bibitte}-genomes/*.genomic.fna.gz \
	-p k=31,scaled=1000,abund \
	--output-dir ${refseq_genomes}/${bibitte}/${dwnld_date}-signatures \
	--singleton # >> singleton means one signature for each record
done

# Create index
$sourmash index -k 31 $refseq_genomes/RefSeq-ABPV-090425.k31.sbt.zip $refseq_genomes/*/${dwnld_date}-signatures/*.sig

cp $refseq_genomes/RefSeq-ABPV-090425.k31.rocksdb

# Extract taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.EXTRA.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.Z

uncompress -c taxdump.tar.Z | tar xf - 
