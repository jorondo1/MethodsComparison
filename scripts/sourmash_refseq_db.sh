ml apptainer

sourmash="singularity exec --writable-tmpfs -e -B $ILAFORES:$ILAFORES,/fast2/def-ilafores:/fast2/def-ilafores $ILL_PIPELINES/containers/sourmash.4.8.11.sif"
refseq_genomes="/fast2/def-ilafores/refseq_genomes"
dwnld_date="20250528"
genomes=$refseq_genomes/${dwnld_date}_genomes
signatures=$refseq_genomes/${dwnld_date}_signatures
mkdir -p $genomes
mkdir -p $signatures

# Extract all accessions + urls from kraken library
cat $ILAFORES/ref_dbs/kraken2_dbs/k2_standard_20241228/library_report.tsv | \
cut -f3 | sed 's|https://ftp.ncbi.nlm.nih.gov/||g' | \
grep genomes | sed 's|genomes/|https://ftp.ncbi.nlm.nih.gov/genomes/|g' | sort -u > $refseq_genomes/ncbi_accession_urls.txt

# DOWNLOAD GENOMES
cd $genomes
xargs -P 8 -I {} bash -c '
    file=$(basename "{}")
    if [ ! -f "$file" ]; then
        curl -L -O "{}" || echo "Failed: {}" >&2
    fi
' < "$refseq_genomes/ncbi_accession_urls.txt"
cd ..

# SKETCH SOURMASH SIGNATURES
# List genomes 
find "$genomes" -name '*.gz' > "$refseq_genomes/${dwnld_date}_genome_list.txt"

# parallel sketch :
nice parallel -j 24 $sourmash sourmash sketch dna {} -p k=31,scaled=1000,abund \
    --name-from-first --outdir "$signatures" \
    :::: "$refseq_genomes/${dwnld_date}_genome_list.txt"

find "$signatures" -name '*.sig' > "$refseq_genomes/${dwnld_date}_signature_list.txt"

# Rename using GCA|F identifier
nice parallel -j 24 'id=$(basename {} | sed '\''s/\.fna\.gz\.sig$//'\''); \
    '"$sourmash"' sourmash signature rename {} "$id" -o {}' \
    :::: "$refseq_genomes/${dwnld_date}_signature_list.txt"

# Find missing signatures 
missing_genomes=($(grep -vFf <(sed "s|$signatures/||g" 20250528_signature_list.txt | sed 's|.sig||g') 20250528_genome_list.txt))
for genome in "${missing_genomes[@]}"; do
    $sourmash sourmash sketch dna "$genome" -p k=31,scaled=1000,abund --name-from-first --outdir "$signatures"
done 
## These genomes seem to be corrupt at the source. Redownloading them did not solve.
###### GCF_045946455.1_ASM4594645v1_genomic.fna.gz
###### GCF_010223795.1_ASM1022379v1_genomic.fna.gz
###### GCF_000847605.1_ViralProj14684_genomic.fna.gz

# BUILD INDEX 
$sourmash sourmash index --dna -k 31 ${signatures}.k31 --from-file $refseq_genomes/${dwnld_date}_signature_list.txt
cp /fast2/def-ilafores/refseq_genomes/RefSeq_20250528.k31.sbt.zip $ILAFORES/ref_dbs/sourmash_db/

# TAXONOMY LINEAGE FILE for sourmash
# 1. Extract just the IDs (GC?_xxxxxxxxx.x
sed -E 's~.*/((GCF|GCA)_[0-9]{9}).*~\1~' $refseq_genomes/ncbi_accession_urls.txt > $refseq_genomes/ncbi_accession.txt

mkdir -p ncbi_ref && cd ncbi_ref

# 2. Download NCBI assembly accession-to-taxid map (updated daily)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt # necessary?

# 3. Download complete taxonomy dump (nodes/names)
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzvf taxdump.tar.gz 
cd ..

# 4. Subset taxids of our database
grep -wFf $refseq_genomes/ncbi_accession.txt <(cat ncbi_ref/assembly_summary_*.txt) | awk -F'\t' '{print $1"\t"$6}' > ncbi_ref/assembly_subset_taxid.tsv

# 5. Create full lineage info
# Requires taxonkit v0.20 since March 2025, NCBI made massive changes to rank names and so on
conda activate
taxonkit reformat2 -I 2 --data-dir ncbi_ref/ \
    -f "{domain};{phylum};{class};{order};{family};{genus};{species}" ncbi_ref/assembly_subset_taxid.tsv \
    | sed 's/\t\|;/,/g' > ncbi_ref/RefSeq_20250528.lineages.csv

# Copy final files to sourmash ref directory
cp ncbi_ref/RefSeq_20250528.lineages.csv $ILAFORES/ref_dbs/sourmash_db/


#find $PWD/split_genomes_09-April-2025/signatures/ -type f -name "*.sig" > split_genomes_09-April-2025/signature_list.txt
#
#singularity exec --writable-tmpfs -e -B $ANCHOR$ILAFORES:$ANCHOR$ILAFORES,$ANCHOR/fast2/def-ilafores:$ANCHOR/fast2/def-ilafores $ANCHOR$ILL_PIPELINES/containers/sourmash.4.8.11.sif sourmash index -k=31 split_genomes_09-April-2025/RefSeq_abpv_09042025 --from-file $ANCHOR$PWD/split_genomes_09-April-2025/signature_list.txt

#
## branchwater multithread sketch
#
## Build csv
#echo name,genome_filename,protein_filename > $refseq_genomes/manysketch_${dwnld_date}.csv
#for i in ${refseq_genomes}/*/${dwnld_date}-*-genomes/*.genomic.fna.gz
#do
#echo $i,$i,
#done >> $refseq_genomes/manysketch_${dwnld_date}.csv
#
#sourmash scripts manysketch fa.csv -o fa.zip -p k=21,k=31,k=51,scaled=1000,abund -p protein,k=10,scaled=200
#
#$sourmash scripts manysketch $refseq_genomes/manysketch_${dwnld_date}.csv \
#--cores 16 -p k=31,scaled=1000,abund --singleton -o $refseq_genomes/all_sig_${dwnld_date}.zip
#
########## OLD SINGLE THREAD
## Compute signature 
#for bibitte in archaea bacteria plasmid viral; do
#	mkdir -p ${refseq_genomes}/${bibitte}/${dwnld_date}-signatures
#	$sourmash sketch dna ${refseq_genomes}/${bibitte}/${dwnld_date}-${bibitte}-genomes/*.genomic.fna.gz \
#	-p k=31,scaled=1000,abund \
#	--output-dir ${refseq_genomes}/${bibitte}/${dwnld_date}-signatures \
#	--singleton # >> singleton means one signature for each record
#done
#
## Create index
#$sourmash sig cat $refseq_genomes/*/${dwnld_date}-signatures/*.sig -o $refseq_genomes/all_sig.zip
## OR index (may take weeks)
#$sourmash index -k 31 $refseq_genomes/RefSeq-ABPV-090425.k31.sbt.zip $refseq_genomes/*/${dwnld_date}-signatures/*.sig
#
#cp $refseq_genomes/RefSeq-ABPV-090425.k31.rocksdb
#
## Extract taxonomy
#wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
#wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
#wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.EXTRA.gz
#wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.Z
#
#uncompress -c taxdump.tar.Z | tar xf - 

