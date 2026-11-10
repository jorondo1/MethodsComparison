#!/bin/bash
# Classification rate per sample, one file per tool in Out/classification_rates/ (no header):
#   path <tab> classified <tab> unclassified <tab> rate   (sourmash, bracken: path <tab> value; reads: path <tab> pairs <tab> reads)
# Source it (pipeline_2_profiling.sh does), then: rate_<tool> [dataset ...]  (no dataset: every data/*/)
# Or run it: bash scripts/src/classification_rates.sh <motus|mpa|kraken|bracken|kraken_species|reads|sourmash|sylph|all> [dataset ...]
# Rows are added or updated by path, never removed: a dataset whose outputs are gone keeps its rows.

# Rows for one dataset, run from $MC (paths start with ./data/<ds>/). awk with ENDFILE = gawk (HPC).
# BEGINFILE, not FNR == 1: an empty profile has no first line and would inherit the previous sample's values
# mOTUs 3 and 4 differ: the count is the last column in both, and the row is unassigned or mOTUv4.0_unassigned
_rows_motus() {
	# -F'\t': taxonomy names contain spaces, so splitting on whitespace puts the columns in the wrong place
	find ./data/"$1"/MOTUS* -name '*_profile.txt' -exec awk -F'\t' '
	BEGINFILE { assigned = 0; unassigned = 0 }
	/^#/ || /^mOTU\t/ { next }
	$1 ~ /unassigned$/ { unassigned = $NF; next }
	{ assigned += $NF }
	ENDFILE {
		total = unassigned + assigned
		if (total > 0) printf "%s\t%.2f\t%.2f\t%.5f\n", FILENAME, assigned, unassigned, assigned / total
		else printf "%s\tNA\tNA\tNA\n", FILENAME
	}' {} + 2> /dev/null
}

# Classified % = sum of the species rows (s__, not t__). With unclassified estimation (--unclassified_estimation before 4.2,
# default since; check the command in line 2 of the profile) abundances are scaled by the mapped fraction, so each rank
# sums to the classified %. The profiles no longer have their UNCLASSIFIED line; 100 - sum(s__) matched it (31.82, AD_Skin SRR20002363).
_rows_mpa() {
	find ./data/"$1"/MPA* -name '*_profile.txt' -exec awk -F'\t' '
	BEGINFILE { classified = 0; lines = 0 }
	!/^#/ { lines++ }
	!/^#/ && $1 ~ /\|s__/ && $1 !~ /\|t__/ { classified += $3 }
	ENDFILE {
		if (lines) printf "%s\t%.2f\t%.2f\t%.5f\n", FILENAME, classified, 100 - classified, classified / 100
		else printf "%s\tNA\tNA\tNA\n", FILENAME
	}' {} + 2> /dev/null
}

# Needs the kraken .kreport files: run before the Kraken cleanup deletes them
_rows_kraken() {
	find ./data/"$1"/KB* -name '*.kreport' ! -name '*bracken*' -exec awk '
	FNR == 1 { unclassified = $2; next }
	FNR == 2 { classified = $2; printf "%s\t%d\t%d\t%.5f\n", FILENAME, classified, unclassified, classified / (unclassified + classified) }
	' {} + 2> /dev/null
}

# Reads at species level after Bracken: path <tab> reads. Divide by the pairs of the reads file (= Kraken's total).
# Kept in the merged tables' source files, so this works after the Kraken cleanup.
_rows_bracken() {
	find ./data/"$1"/KB* -name '*_bracken_S.MPA.TXT' -exec awk -F'\t' '
	BEGINFILE { reads = 0 }
	!/^#/ && $1 ~ /\|s__/ { reads += $NF }
	ENDFILE { printf "%s\t%d\n", FILENAME, reads }
	' {} + 2> /dev/null
}

# Species-level identifiability of Kraken: reads Kraken itself assigned at species rank or below, over all reads.
# Not the Bracken total (bracken_S.MPA.TXT), which redistributes nearly every classified read down to species, from any
# rank. Read from the BRACKEN SUMMARY that kraken_local.sh tees into <KB dir>/kraken_wrapper.log, so no .kreport needed:
#   species = "kept at species level" + "discarded (species reads < threshold)"; total = "Total reads in sample".
# Rows: ./data/<ds>/<KB>/<sample>/<sample>.kreport (same path layout as the kraken file) <tab> species <tab> other <tab> rate.
# Only samples of the dataset's sample tsv; a log can hold several runs of a sample: the last summary wins.
_rows_kraken_species() {
	local tsv log kb
	tsv=$(awk -F'\t' -v d="$1" '$1 == d { print $2 }' ./scripts/src/datasets.tsv)
	tsv=$(eval echo "$tsv")
	[[ -f $tsv ]] || return 0
	for log in ./data/"$1"/KB*/kraken_wrapper.log; do
		[[ -f $log ]] || continue
		kb=$(basename "$(dirname "$log")")
		awk -v pre="./data/$1/$kb" '
		NR == FNR { keep[$1]; next }
		/BRACKEN SUMMARY \(Kraken report:/ { s = $0; sub(/\.kreport\).*/, "", s); sub(/.*\//, "", s); next }
		/Total reads in sample:/ { tot[s] = $NF + 0; sp[s] = 0; next }
		/kept at species level/ || /discarded \(species reads/ { sp[s] += $NF; next }
		END {
			for (x in tot) {
				if (!(x in keep)) continue
				if (tot[x] > 0) printf "%s/%s/%s.kreport\t%d\t%d\t%.5f\n", pre, x, x, sp[x], tot[x] - sp[x], sp[x] / tot[x]
				else printf "%s/%s/%s.kreport\tNA\tNA\tNA\n", pre, x, x
			}
		}' "$tsv" "$log"
	done
}

# Reads per sample after KneadData (final counts in its log): ./data/<ds>/preproc/<sample> <tab> pairs <tab> reads.
# pairs = Kraken's total (only the paired reads go to Kraken), so it is the Bracken denominator even where the kreports
# are gone. reads = pairs x 2 + orphans: the depth shown in the identifiability figure. MetaPhlAn's "reads processed"
# is ~1 % lower (it drops reads under 70 bp). A log can hold several runs: the last count wins.
_rows_reads() {
	local tsv sample p1 log
	tsv=$(awk -F'\t' -v d="$1" '$1 == d { print $2 }' ./scripts/src/datasets.tsv)
	tsv=$(eval echo "$tsv")
	[[ -f $tsv ]] || return 0
	while IFS=$'\t' read -r sample p1 _; do
		log=$(dirname "$p1")/kneaddata-$sample.log
		# sample tsv written on another cluster (Horse: Fir paths): the log copied into this repo's data/<ds>/preproc/<sample>/
		[[ -f $log ]] || log=./data/$1/preproc/$sample/kneaddata-$sample.log
		[[ -f $log ]] || continue
		awk -v key="./data/$1/preproc/$sample" '
		match($0, /READ COUNT: final (pair1|orphan1|orphan2) :.*: ([0-9.]+)$/, m) { n[m[1]] = m[2] + 0 }
		END { if ("pair1" in n) printf "%s\t%d\t%d\n", key, n["pair1"], 2 * n["pair1"] + n["orphan1"] + n["orphan2"] }
		' "$log"
	done < "$tsv"
}

_rows_sourmash() {
	find ./data/"$1"/SM_* -name '*_gather.csv' -exec awk -F',' '
	FNR == 1 {
		for (i = 1; i <= NF; i++) if ($i == "f_unique_weighted") { col = i; break }
		next
	}
	{ sum += $col }
	END { if (col) printf "%s\t%.5f\n", FILENAME, sum }
	' {} \; 2> /dev/null
}

# Written by sylph_merge.sh
_rows_sylph() {
	cat ./data/"$1"/SYLPH_*/classification_rate.tsv 2> /dev/null
}

# _update_rates <tool> [dataset ...]
_update_rates() {
	local tool=$1; shift
	local mc=${MC:-$PWD}
	local file=$mc/Out/classification_rates/${tool}_classification_rate.tsv
	local ds rows tmp n_old n_kept
	[[ -d $mc/data ]] || { echo "No data/ in $mc: set \$MC or run from the repo"; return 1; }
	(( $# )) || set -- $(cd "$mc" && ls -d data/*/ | xargs -n1 basename)
	mkdir -p "$(dirname "$file")" && touch "$file"
	for ds in "$@"; do
		rows=$(cd "$mc" && "_rows_$tool" "$ds")
		n_old=$(grep -c "^\./data/$ds/" "$file")
		if [[ -z $rows ]]; then
			echo "$tool $ds: no output files, $n_old existing rows kept"
			continue
		fi
		tmp=$(mktemp)
		# keep the old rows whose path isn't in the new rows, then add the new rows
		awk -F'\t' 'NR == FNR { new[$1]; next } !($1 in new)' <(echo "$rows") "$file" > "$tmp"
		n_kept=$(grep -c "^\./data/$ds/" "$tmp")
		echo "$rows" >> "$tmp"
		sort "$tmp" > "$file" && rm "$tmp"
		echo "$tool $ds: $(wc -l <<< "$rows" | tr -d ' ') rows written, $n_kept old rows without a new value kept"
	done
}

rate_motus() { _update_rates motus "$@"; }
rate_mpa() { _update_rates mpa "$@"; }
rate_kraken() { _update_rates kraken "$@"; }
rate_bracken() { _update_rates bracken "$@"; }
rate_kraken_species() { _update_rates kraken_species "$@"; }
rate_reads() { _update_rates reads "$@"; }
rate_sourmash() { _update_rates sourmash "$@"; }
rate_sylph() { _update_rates sylph "$@"; }
rate_all() { local t; for t in motus mpa kraken bracken kraken_species reads sourmash sylph; do _update_rates $t "$@"; done; }

if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
	tool=${1:?usage: bash classification_rates.sh <motus|mpa|kraken|bracken|kraken_species|reads|sourmash|sylph|all> [dataset ...]}; shift
	"rate_$tool" "$@"
fi
