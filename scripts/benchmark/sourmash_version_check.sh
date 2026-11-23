#!/bin/bash
# Are the sourmash benchmark profiles the same with 4.9.4 + branchwater 0.9.15 (fastgather, RefSeq zip) as with
# 4.8.11 (RefSeq: Python gather on the SBT; GTDB: branchwater 0.9.5 fastgather)? Inputs and tax step, part 1 of 2.
#   1. sketch md5 of every benchmark sample, both versions (sig describe)
#   2. species summaries of the four gather sets with the production tax step: sourmash_tax_summary.py in the 4.8.11
#      container for all of them (same tax code), same lineages as merge_sourmash() in merge_tables.sh
# Then: Rscript scripts/benchmark/compare_sourmash_versions.R (metrics, verdict).
# Usage, from the repo root: bash scripts/benchmark/sourmash_version_check.sh   (light: OK on the interactive node)
set -eo pipefail

MC=$ILAFORES/analysis/MethodsComparison
cd $MC
B=Out/benchmark
OUT=$B/version_check
SM_DB=$ILAFORES/ref_dbs/sourmash_db
SYLPH_REFSEQ_TAX=$ILAFORES/ref_dbs/sylph_db/RefSeq_20250528_taxonomy.tsv.gz
CONTAINERS=${ILL_PIPELINES:-$ILAFORES/programs/ILL_pipelines}/containers
mkdir -p $OUT
module load StdEnv/2020 apptainer/1.1.5
SING="singularity exec -e -B $ILAFORES"

# 1. Sketches. The 4.8.11 sketches in signatures/, the 4.9.4 ones in sm49/signatures/
for v in 4.8.11:$B/signatures 4.9.4:$B/sm49/signatures; do
	ver=${v%%:*} dir=${v#*:}
	$SING $CONTAINERS/sourmash.$ver.sif sourmash sig describe --csv $OUT/sigs_$ver.csv $dir/*.sig > /dev/null 2>&1
	echo "sketches $ver: $(($(wc -l < $OUT/sigs_$ver.csv) - 1))"
done

# 2. RefSeq lineages: same awk as merge_sourmash() (ident = genome file name without extension = match name)
zcat "$SYLPH_REFSEQ_TAX" | awk -F'\t' '
	BEGIN { print "ident,superkingdom,phylum,class,order,family,genus,species" }
	{
		id = $1; sub(/\.(fna|fa|fasta)(\.gz)?$/, "", id)
		n = split($2, r, ";"); line = id
		for (i = 1; i <= 7; i++) { x = r[i]; sub(/^[a-z]__/, "", x); if (x ~ /[,"]/) { gsub(/"/, "\"\"", x); x = "\"" x "\"" }; line = line "," x }
		print line
	}' > $OUT/refseq.lineages.csv

# 3. Tax summaries. Gather CSVs with a header only (no match above 50 kbp) can't be loaded: listed, left out
for set in old:$B/SM_RefSeq_20250528 new:$B/sm49/SM_RefSeq_20250528 old:$B/SM_gtdb-rs220-reps new:$B/sm49/SM_gtdb-rs220-reps; do
	ver=${set%%:*} dir=${set#*:} db=${dir##*/SM_}
	case $db in
		RefSeq_20250528) lin=($OUT/refseq.lineages.csv); opts=(--keep-full-identifiers --keep-identifier-versions) ;;
		gtdb-rs220-reps) lin=($SM_DB/gtdb-rs220.lineages.reps.csv); opts=() ;;
	esac
	list=$OUT/${db}_$ver.list
	: > $list
	for f in $dir/*_gather.csv; do
		if (( $(wc -l < $f) > 1 )); then echo $MC/$f >> $list; else echo "  $db $ver: no match in $(basename $f)"; fi
	done
	$SING $CONTAINERS/sourmash.4.8.11.sif python3 -I $MC/scripts/src/sourmash_tax_summary.py --from-file $list \
		-t "${lin[@]}" "${opts[@]}" -o $OUT/${db}_$ver.summary.csv > $OUT/${db}_$ver.log 2>&1 \
		|| { tail -5 $OUT/${db}_$ver.log; exit 1; }
	echo "$db $ver: $(wc -l < $list) gather files, $(grep -c ',species,' $OUT/${db}_$ver.summary.csv) species rows;" \
		"$(grep -o 'lineage assignments for [0-9]* results were missed' $OUT/${db}_$ver.log | awk '{s += $4} END {print s + 0}') matches without lineage"
done
gzip -fn $OUT/*.summary.csv
