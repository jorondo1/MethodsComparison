# sourmash RefSeq database (2025-05-28): the genomes of the Kraken2 standard 20241228 library, sketched with sourmash,
# indexed, and given NCBI lineages. Log of the commands used.
# Copy-paste collection: running or sourcing this file does nothing. Paste a block, wait for it to finish, then the next.
# Older approaches (manysketch, per-kingdom singleton sketches, accession2taxid): git history of this file.
return 0 2>/dev/null; [[ $0 == "${BASH_SOURCE[0]}" ]] && { echo "Copy-paste collection: paste the commands, don't run this file"; exit 1; }

module load StdEnv/2020 apptainer/1.1.5
sourmash="singularity exec --writable-tmpfs -e -B $ILAFORES:$ILAFORES,/fast2/def-ilafores:/fast2/def-ilafores $ILL_PIPELINES/containers/sourmash.4.8.11.sif"
refseq_genomes=/fast2/def-ilafores/refseq_genomes
dwnld_date=20250528
genomes=$refseq_genomes/${dwnld_date}_genomes
signatures=$refseq_genomes/${dwnld_date}_signatures
SM_DB=$ILAFORES/ref_dbs/sourmash_db
mkdir -p $genomes $signatures
cd $refseq_genomes

################################################################
# Genomes ########
################################################################

# Accessions + URLs of the Kraken2 standard library
cut -f3 $ILAFORES/ref_dbs/kraken2_dbs/k2_standard_20241228/library_report.tsv | \
	sed 's|https://ftp.ncbi.nlm.nih.gov/||g' | grep genomes | \
	sed 's|genomes/|https://ftp.ncbi.nlm.nih.gov/genomes/|g' | sort -u > ncbi_accession_urls.txt

# Download (skips files already there)
cd $genomes
xargs -P 8 -I {} bash -c '
	file=$(basename "{}")
	[ -f "$file" ] || curl -L -O "{}" || echo "Failed: {}" >&2
' < $refseq_genomes/ncbi_accession_urls.txt
cd $refseq_genomes
find $genomes -name '*.gz' > ${dwnld_date}_genome_list.txt

################################################################
# Signatures: one per genome, k=31, scaled=1000, with abundances ########
################################################################

nice parallel -j 24 $sourmash sourmash sketch dna {} -p k=31,scaled=1000,abund \
	--name-from-first --outdir $signatures :::: ${dwnld_date}_genome_list.txt
find $signatures -name '*.sig' > ${dwnld_date}_signature_list.txt

# Name each signature after its genome file (GCF_/GCA_ accession + assembly name): this is the gather match name,
# which merge_tables.sh joins to the lineages
nice parallel -j 24 'id=$(basename {} | sed '\''s/\.fna\.gz\.sig$//'\''); \
	'"$sourmash"' sourmash signature rename {} "$id" -o {}' :::: ${dwnld_date}_signature_list.txt

# Sketch the genomes that have no signature yet
missing_genomes=($(grep -vFf <(sed "s|$signatures/||; s|\.sig$||" ${dwnld_date}_signature_list.txt) ${dwnld_date}_genome_list.txt))
for genome in "${missing_genomes[@]}"; do
	$sourmash sourmash sketch dna "$genome" -p k=31,scaled=1000,abund --name-from-first --outdir $signatures
done
# These 3 stay missing: corrupt at the source, and downloading them again didn't help.
#   GCF_045946455.1_ASM4594645v1_genomic.fna.gz
#   GCF_010223795.1_ASM1022379v1_genomic.fna.gz
#   GCF_000847605.1_ViralProj14684_genomic.fna.gz
find $signatures -name '*.sig' > ${dwnld_date}_signature_list.txt

################################################################
# Index ########
################################################################

# SBT index (2025-06), used by the main analysis (sourmash 4.8.11 gather, 1 thread)
$sourmash sourmash index --dna -k 31 RefSeq_${dwnld_date}.k31 --from-file ${dwnld_date}_signature_list.txt
cp RefSeq_${dwnld_date}.k31.sbt.zip $SM_DB/

# Zip collection (2026-11-19), for branchwater fastgather (multithreaded), which can't read an SBT. Same signatures:
# sig cat copies the SBT's leaves. The output name must end in .zip (sourmash picks the format from the extension).
sbatch -A def-ilafores -p bio -J sm_refseq_zip -c 1 --mem=64G --time=12:00:00 -o sm_refseq_zip-%j.out --wrap \
	"module load StdEnv/2020 apptainer/1.1.5; $sourmash sourmash sig cat $SM_DB/RefSeq_${dwnld_date}.k31.sbt.zip \
	-o $SM_DB/RefSeq_${dwnld_date}.k31.tmp.zip && mv $SM_DB/RefSeq_${dwnld_date}.k31.tmp.zip $SM_DB/RefSeq_${dwnld_date}.k31.zip"
# Check: same number of sketches, ksize, scaled, abundances and total hashes in both
$sourmash sourmash sig fileinfo $SM_DB/RefSeq_${dwnld_date}.k31.sbt.zip
$sourmash sourmash sig fileinfo $SM_DB/RefSeq_${dwnld_date}.k31.zip

################################################################
# Lineages ########
################################################################

# Note: merge_tables.sh no longer uses this file for sourmash. It lists GCA/GCF twins and several versions of some
# accessions, with no header, so it uses the sylph-tax RefSeq table instead (merge_sourmash() in merge_tables.sh).

# 1. Accessions (GCF_/GCA_ + 9 digits)
sed -E 's~.*/((GCF|GCA)_[0-9]{9}).*~\1~' ncbi_accession_urls.txt > ncbi_accession.txt

# 2. NCBI assembly summaries and taxonomy dump. The historical summaries list suppressed assemblies, needed to
#    reproduce the 20241228 Kraken library.
mkdir -p ncbi_ref && cd ncbi_ref
for f in refseq genbank genbank_historical refseq_historical; do
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_$f.txt
done
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzvf taxdump.tar.gz
cd $refseq_genomes

# 3. Taxid of each accession of the database
grep -wFf ncbi_accession.txt <(cat ncbi_ref/assembly_summary_*.txt) | awk -F'\t' '{print $1"\t"$6}' > ncbi_ref/assembly_subset_taxid.tsv
# Accessions whose GCF_ was later removed from RefSeq are only in GenBank, as GCA_: match them on the number.
# (Probably unnecessary now that the historical summaries are included.)
cut -f1 ncbi_ref/assembly_subset_taxid.tsv | sed 's/^[^_]*_\([^.]*\).*$/\1/' > ncbi_ref/acc_nums.txt
missing_num_acc=($(grep -vFf ncbi_ref/acc_nums.txt ncbi_accession.txt | grep GC | sed 's|GC.||g'))
grep -f <(printf "%s\n" "${missing_num_acc[@]}") <(cat ncbi_ref/assembly_summary_*.txt) | \
	awk -F'\t' '{print $1"\t"$6}' >> ncbi_ref/assembly_subset_taxid.tsv
# A few remain without a taxid (e.g. GCA_045946455, suppressed by NCBI).

# 4. Lineages. taxonkit >= 0.20: NCBI renamed ranks in March 2025.
conda activate
taxonkit reformat2 -I 2 --data-dir ncbi_ref/ \
	-f "{domain};{phylum};{class};{order};{family};{genus};{species}" ncbi_ref/assembly_subset_taxid.tsv \
	| sed 's/\t\|;/,/g' > ncbi_ref/RefSeq_${dwnld_date}.lineages.csv
cp ncbi_ref/RefSeq_${dwnld_date}.lineages.csv $SM_DB/
