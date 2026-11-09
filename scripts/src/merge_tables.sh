#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH -D /net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison
#SBATCH -o /net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison/logs/merge_tables-%j.slurm.out
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -A def-ilafores
#SBATCH -J merge_tables

# Merge each profiler's per-sample outputs into one table per dataset x database, with the tool's own merge utility.
# R (scripts/0_proc_data.R) reads these tables instead of the per-sample files. Per-sample files are never modified.
#
# Usage (from $MC, or sbatch): merge_tables.sh [-f] [--tools motus,mpa,kraken,sourmash,sylph] [dataset ...]
#   no dataset   : every dataset of scripts/src/datasets.tsv that has a data/<dataset>/ folder
#   --tools      : comma-separated subset of motus,mpa,kraken,sourmash,sylph (default: all)
#   -f, --force  : merge a folder even when samples of the dataset's sample tsv have no output file (they're listed);
#                  without -f such folders are skipped. Passed on to sylph_merge.sh.
#
# Outputs, next to the per-sample files:
#   data/<ds>/MOTUS/merged_profile.tsv          mOTUs 3.1.0: own outer join (motus merge 3.1.0 needs identical rows in
#                                               every profile, but motus3_SLURM.sh drops the zero rows). Empty profiles skipped.
#   data/<ds>/MOTUS4/merged_profile.tsv         mOTUs 4.0.4: motus merge
#   data/<ds>/MPA_<db>/merged_profile.tsv       MetaPhlAn 4.2.4: merge_metaphlan_tables.py
#   data/<ds>/KB<x>/merged_profile.tsv          KrakenTools 1.2: combine_mpa.py on the *_bracken_S.MPA.TXT files
#                                               (header rewritten with the sample names: the inputs have no header line)
#   data/<ds>/SM_<db>/merged_tax_summary.csv.gz sourmash 4.8.11: tax metagenome csv_summary, all samples
#                                               (sourmash_tax_summary.py: the CLI's csv_summary takes one sample at a time)
#   data/<ds>/SYLPH_<db>/...                    sylph_merge.sh (counts.tsv, taxonomy.tsv, merged_species.tsv, classification_rate.tsv)
# Sample name = per-sample file name minus its suffix, as in community_functions.R read_filename()/parse_SM().
set -eo pipefail

FORCE=0
TOOLS=motus,mpa,kraken,sourmash,sylph
DATASETS=()
while [[ $# -gt 0 ]]; do
	case $1 in
		-f|--force) FORCE=1 ;;
		--tools) TOOLS=$2; shift ;;
		-h|--help) sed -n '/^# Merge each/,/^set -eo/p' "$0" | sed '$d'; exit 0 ;;
		-*) echo "Unknown option $1"; exit 1 ;;
		*) DATASETS+=("$1") ;;
	esac
	shift
done

MC=${MC:-$PWD}
[[ -d $MC/data && -f $MC/scripts/src/datasets.tsv ]] || { echo "Run from the MethodsComparison repo (or set \$MC)"; exit 1; }
cd "$MC"

ILAFORES=${ILAFORES:-/net/nfs-bio/jbod2/def-ilafores}
CONTAINERS=${ILL_PIPELINES:-$ILAFORES/programs/ILL_pipelines}/containers
MPA_SIF=$CONTAINERS/metaphlan.4.2.4.sif
MOTUS4_SIF=$CONTAINERS/mOTUs_v4.0.4.sif
SM_SIF=$CONTAINERS/sourmash.4.8.11.sif
COMBINE_MPA=$ILAFORES/programs/KrakenTools-1.2/combine_mpa.py
SM_DB=$ILAFORES/ref_dbs/sourmash_db
SYLPH_REFSEQ_TAX=$ILAFORES/ref_dbs/sylph_db/RefSeq_20250528_taxonomy.tsv.gz

command -v singularity > /dev/null || module load apptainer 2> /dev/null || true
command -v singularity > /dev/null || { echo "singularity not found on $(hostname)"; exit 1; }
WORK=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/merge_tables.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

SING="singularity exec -e -B $ILAFORES -B $WORK" # $WORK: lists and lineage files passed to the containers
[[ $MC == "$ILAFORES"/* ]] || SING="$SING -B $MC"


has_tool() { [[ ",$TOOLS," == *",$1,"* ]]; }

declare -A SAMPLE_TSV
while IFS=$'\t' read -r ds tsv; do
	[[ -z $ds || $ds == \#* ]] && continue
	SAMPLE_TSV[$ds]=$(eval echo "$tsv")
done < scripts/src/datasets.tsv

if [[ ${#DATASETS[@]} -eq 0 ]]; then
	while IFS=$'\t' read -r ds tsv; do
		[[ -z $ds || $ds == \#* ]] && continue
		[[ -d data/$ds ]] && DATASETS+=("$ds")
	done < scripts/src/datasets.tsv
fi

FAILED=()

# check_samples <suffix> <files...>: samples of the tsv without an output file. Returns 1 if the folder must be skipped.
check_samples() {
	local suffix=$1; shift
	check_names "$(for f in "$@"; do basename "$f" "$suffix"; done)"
}

# check_names <newline-separated sample names>
check_names() {
	local tsv=${SAMPLE_TSV[$DS]}
	[[ -f $tsv ]] || { echo "  no sample tsv for $DS, can't check for missing samples"; return 0; }
	local found missing extra
	found=$(echo "$1" | sort -u)
	missing=$(comm -23 <(cut -f1 "$tsv" | sort -u) <(echo "$found"))
	extra=$(comm -13 <(cut -f1 "$tsv" | sort -u) <(echo "$found") | wc -l)
	[[ $extra -gt 0 ]] && echo "  $extra samples not in the sample tsv (merged anyway)"
	if [[ -n $missing ]]; then
		echo "  $(echo "$missing" | wc -l) samples of the tsv have no output: $(echo $missing)"
		if [[ $FORCE -eq 0 ]]; then echo "  skipped (-f to merge without them)"; return 1; fi
	fi
	return 0
}

# check_columns <merged file> <header pattern> <n leading columns> <suffix> <files...>: sample columns = file names
# (as a set: motus merge sorts the columns)
check_columns() {
	local merged=$1 pattern=$2 lead=$3 suffix=$4; shift 4
	diff -q <(grep -m1 "$pattern" "$merged" | cut -f$((lead + 1))- | tr '\t' '\n' | sort) \
		<(for f in "$@"; do basename "$f" "$suffix"; done | sort) > /dev/null \
		|| { echo "  sample columns of $merged don't match the input file names"; return 1; }
}

# finish <tmp file> <final file>: report and move into place
finish() {
	mv "$1" "$2"
	echo "  -> $2"
}

merge_motus3() {
	local dir=$1 out=$1/merged_profile.tsv files=() empty=()
	for f in "$dir"/*_profile.txt; do [[ -s $f ]] && files+=("$f") || empty+=("$(basename "$f" _profile.txt)"); done
	[[ ${#empty[@]} -gt 0 ]] && echo "  ${#empty[@]} empty profiles skipped (also absent from the per-sample parsing): ${empty[*]}"
	[[ ${#files[@]} -gt 0 ]] || { echo "  no profile"; return 1; }
	check_samples _profile.txt "${files[@]}" "${empty[@]/%/_profile.txt}" || return 1
	awk -F'\t' -v OFS='\t' '
		FNR == 1 { n++; s = FILENAME; sub(/.*\//, "", s); sub(/_profile\.txt$/, "", s); name[n] = s; delete here }
		/^#/ { if (n == 1 && FNR == 1) version = $0; next }
		{
			key = $1 OFS $2 OFS $3
			if (key in here) { print "duplicate row in " FILENAME ": " key > "/dev/stderr"; exit 1 }
			here[key] = 1
			if (!(key in row)) { row[key] = ++k; keys[k] = key }
			val[row[key], n] = $4
		}
		END {
			print version
			print "# merged by merge_tables.sh: outer join of " n " mOTUs 3 profiles (rows missing from a profile = 0)"
			printf "mOTU\tconsensus_taxonomy\tNCBI_tax_id"; for (j = 1; j <= n; j++) printf "\t%s", name[j]; printf "\n"
			for (i = 1; i <= k; i++) {
				printf "%s", keys[i]
				for (j = 1; j <= n; j++) printf "\t%s", ((i, j) in val) ? val[i, j] : 0
				printf "\n"
			}
		}' "${files[@]}" > "$WORK/m3.tsv"
	check_columns "$WORK/m3.tsv" '^mOTU' 3 _profile.txt "${files[@]}" || return 1
	echo "  ${#files[@]} samples, $(awk '/^mOTU\t/ {h = 1; next} h' "$WORK/m3.tsv" | wc -l) mOTUs"
	finish "$WORK/m3.tsv" "$out"
}

merge_motus4() {
	local dir=$1 out=$1/merged_profile.tsv files=("$1"/*_profile.txt)
	[[ -f ${files[0]} ]] || { echo "  no profile"; return 1; }
	check_samples _profile.txt "${files[@]}" || return 1
	printf '%s\n' "${files[@]/#/$MC/}" > "$WORK/m4.list"
	$SING "$MOTUS4_SIF" motus merge -i "$WORK/m4.list" -o "$WORK/m4.tsv" > "$WORK/m4.log" 2>&1 \
		|| { tail -5 "$WORK/m4.log"; return 1; }
	check_columns "$WORK/m4.tsv" '^mOTU' 2 _profile.txt "${files[@]}" || return 1
	echo "  ${#files[@]} samples, $(awk '/^mOTU\t/ {h = 1; next} h' "$WORK/m4.tsv" | wc -l) mOTUs"
	finish "$WORK/m4.tsv" "$out"
}

merge_mpa() {
	local dir=$1 out=$1/merged_profile.tsv files=("$1"/*/*_profile.txt)
	[[ -f ${files[0]} ]] || { echo "  no profile"; return 1; }
	check_samples _profile.txt "${files[@]}" || return 1
	printf '%s\n' "${files[@]/#/$MC/}" > "$WORK/mpa.list"
	rm -f "$WORK/mpa.tsv"
	$SING "$MPA_SIF" merge_metaphlan_tables.py -l "$WORK/mpa.list" -o "$WORK/mpa.tsv" > "$WORK/mpa.log" 2>&1 \
		|| { tail -5 "$WORK/mpa.log"; return 1; }
	check_columns "$WORK/mpa.tsv" '^clade_name' 1 _profile.txt "${files[@]}" || return 1
	echo "  ${#files[@]} samples, $(grep -vc '^#\|^clade_name' "$WORK/mpa.tsv") clades"
	finish "$WORK/mpa.tsv" "$out"
}

merge_kraken() {
	local dir=$1 out=$1/merged_profile.tsv files=("$1"/*/*_bracken/*_bracken_S.MPA.TXT)
	[[ -f ${files[0]} ]] || { echo "  no *_bracken_S.MPA.TXT"; return 1; }
	check_samples _bracken_S.MPA.TXT "${files[@]}" || return 1
	python3 "$COMBINE_MPA" -i "${files[@]}" -o "$WORK/kb.raw" > "$WORK/kb.log" 2>&1 || { tail -5 "$WORK/kb.log"; return 1; }
	# combine_mpa.py names the columns "Sample #<i>" (input order) when the inputs have no header line
	{ printf '#Classification'; for f in "${files[@]}"; do printf '\t%s' "$(basename "$f" _bracken_S.MPA.TXT)"; done; printf '\n'
	  # combine_mpa.py prints a taxon twice (same values) when it's a root in one input and, in another, comes after a line
	  # that is a prefix of it (inputs with higher ranks, e.g. some RA_Gut KB10 files): keep the first copy
	  tail -n +2 "$WORK/kb.raw" | awk '!seen[$0]++'; } > "$WORK/kb.tsv"
	local dup; dup=$(tail -n +2 "$WORK/kb.tsv" | cut -f1 | sort | uniq -d | wc -l)
	[[ $dup -eq 0 ]] || { echo "  $dup taxa appear twice with different values"; return 1; }
	echo "  ${#files[@]} samples, $(tail -n +2 "$WORK/kb.tsv" | wc -l) taxa"
	finish "$WORK/kb.tsv" "$out"
}

merge_sourmash() {
	local dir=$1 db=${1##*/SM_} out=$1/merged_tax_summary.csv.gz files=("$1"/*_gather.csv) lin=() opts=()
	case $db in
		genbank-2022.03) lin=("$SM_DB"/genbank-2022.03-*.lineages.csv) ;;
		gtdb-rs214-full|gtdb-rs214-rep) lin=("$SM_DB/gtdb-rs214.lineages.csv") ;;
		gtdb-rs220-rep|gtdb-rs220-reps) lin=("$SM_DB/gtdb-rs220.lineages.reps.csv") ;;
		RefSeq_20250528)
			# RefSeq_20250528.lineages.csv has no header and lists GCA/GCF twins and several versions of an accession,
			# some with different lineages, so sourmash can't read it (and R's join on the accession number counted most
			# genomes 2+ times). Use the sylph-tax table instead: one row per genome file, built for the same genome list.
			# ident = genome file name without extension = the gather match name (--keep-full-identifiers).
			[[ -f $WORK/refseq.lineages.csv ]] || zcat "$SYLPH_REFSEQ_TAX" | awk -F'\t' '
				BEGIN { print "ident,superkingdom,phylum,class,order,family,genus,species" }
				{
					id = $1; sub(/\.(fna|fa|fasta)(\.gz)?$/, "", id)
					n = split($2, r, ";"); line = id
					for (i = 1; i <= 7; i++) { x = r[i]; sub(/^[a-z]__/, "", x); if (x ~ /[,"]/) { gsub(/"/, "\"\"", x); x = "\"" x "\"" }; line = line "," x }
					print line
				}' > "$WORK/refseq.lineages.csv"
			lin=("$WORK/refseq.lineages.csv"); opts=(--keep-full-identifiers --keep-identifier-versions) ;;
		gtdb_rs220) echo "  full GTDB rs220 gather: not used, skipped"; return 0 ;;
		*) echo "  no lineage file known for $db, skipped"; return 1 ;;
	esac
	[[ -f ${files[0]} ]] || { echo "  no *_gather.csv"; return 1; }
	# Truncated gather CSVs (an interrupted write, e.g. AD_Skin SRR20002431 rs214-full) end with an incomplete row that
	# sourmash can't load. Those rows are dropped from a copy in $WORK; parse_SM() put such a row under a bogus sample
	# name (the cut query_name), which assemble_phyloseq() dropped, so the result is the same.
	python3 -I - "$WORK/sm_clean" "${files[@]/#/$MC/}" > "$WORK/sm.list" <<'EOF'
import csv, os, sys
out = sys.argv[1]
os.makedirs(out, exist_ok=True)
for f in sys.argv[2:]:
    with open(f, newline='') as fh:
        rows = list(csv.reader(fh))
    bad = [i for i, r in enumerate(rows[1:], 2) if len(r) != len(rows[0])]
    if bad:
        p = os.path.join(out, os.path.basename(f))
        with open(p, 'w', newline='') as o:
            csv.writer(o, lineterminator='\n').writerows(r for r in rows if len(r) == len(rows[0]))
        print(f"  {os.path.basename(f)}: incomplete row(s) {bad} dropped (truncated file?)", file=sys.stderr)
        print(p)
    else:
        print(f)
EOF
	$SING "$SM_SIF" python3 -I "$MC/scripts/src/sourmash_tax_summary.py" --from-file "$WORK/sm.list" \
		-t "${lin[@]}" "${opts[@]}" -o "$WORK/sm.csv" > "$WORK/sm.log" 2>&1 || { tail -5 "$WORK/sm.log"; return 1; }
	local missed; missed=$(grep -o 'lineage assignments for [0-9]* results were missed' "$WORK/sm.log" | awk '{s += $4} END {print s + 0}')
	[[ $missed -gt 0 ]] && echo "  $missed gather matches without lineage (reported as unclassified): $(grep -o 'missing from the taxonomy information: .*' "$WORK/sm.log" | cut -d: -f2 | tr ',' '\n' | sort -u | wc -l) distinct"
	# sample = query_name (what parse_SM() uses); gather file names don't follow one pattern
	local queries; queries=$(tail -n +2 "$WORK/sm.csv" | cut -d, -f1 | sort -u)
	echo "  ${#files[@]} gather files, $(echo "$queries" | wc -l) queries, $(grep -c ',species,' "$WORK/sm.csv") species rows"
	check_names "$queries" || return 1
	gzip -n "$WORK/sm.csv"
	finish "$WORK/sm.csv.gz" "$out"
}

merge_sylph() {
	local dir=$1 db=${1##*/SYLPH_} f=()
	[[ -d $dir/per_sample ]] || { echo "  no per_sample/"; return 1; }
	[[ $FORCE -eq 1 ]] && f=(-f)
	bash scripts/utils/sylph_merge.sh "${f[@]}" "$DS" "$db" "${SAMPLE_TSV[$DS]}" 2>&1 | sed 's/^/  /'
	return "${PIPESTATUS[0]}"
}

run() { # run <tool fn> <dir>
	echo "== $DS/$(basename "$2")"
	"$1" "$2" || FAILED+=("$DS/$(basename "$2")")
}

for DS in "${DATASETS[@]}"; do
	[[ -d data/$DS ]] || { echo "data/$DS not found"; FAILED+=("$DS"); continue; }
	for d in data/"$DS"/*/; do
		d=${d%/}; name=$(basename "$d")
		case $name in
			MOTUS)   has_tool motus    && run merge_motus3   "$d" ;;
			MOTUS4)  has_tool motus    && run merge_motus4   "$d" ;;
			MPA_*)   has_tool mpa      && run merge_mpa      "$d" ;;
			KB*)     has_tool kraken   && run merge_kraken   "$d" ;;
			SM_*)    has_tool sourmash && run merge_sourmash "$d" ;;
			SYLPH_*) has_tool sylph    && run merge_sylph    "$d" ;;
		esac
	done
done

if [[ ${#FAILED[@]} -gt 0 ]]; then
	echo; echo "Not merged (see above): ${FAILED[*]}"; exit 1
fi
echo; echo "All merged."
