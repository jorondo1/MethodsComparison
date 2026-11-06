#!/bin/bash

# Verify output for all datasets
# Sets $FOUND: the tsv rows (array indices) of the samples missing for the last database checked
check_output() {
	if [ $# -lt 3 ]; then
	        echo "!!Usage: check_output <'database1 database2 ... databaseN'> <dataset path> <filename_suffix>"
	        echo "Depends on variables set by use_dataset."
	        return 1  # Exit the function with a non-zero status
	fi

	DATABASES="$1"
	DATASET_PATH="$2"
	filename_suffix="$3"

    if [[ -z "$DATASET_PATH" || -z "$TSV" ]]; then
		echo "\$DATASET_PATH or \$TSV not set. Run use_dataset <dataset> first."
		return 1
	fi

	for db in $DATABASES; do
		FOUND=""
		found=$(find $DATASET_PATH/*$db -type f -name "*$filename_suffix" -exec basename {} \; | sed "s/${filename_suffix}//" | sed "s/_${db}//")

		num_found=$(echo $found | wc -w)
		echo "$num_found $db output for $DATASET found, $N_SAMPLES expected."

		if [ "$num_found" -lt "$N_SAMPLES" ]; then
			export FOUND=$(grep -nv -E "${found}" "$TSV" | cut -d: -f1| paste -s -d,)
			echo $FOUND
		fi
	done
}

# check_all [-m] [dataset ...]: table of ok/empty/unfinished samples for every tool x database x dataset
# (all datasets in scripts/src/datasets.tsv if none given). Samples come from column 1 of each sample tsv.
#   ok = output with results; empty = output with no results; unfinished = no output for that sample
#   Red: unfinished > 0. Yellow: only empty. Green: all ok. No colour when the output isn't a terminal.
#   -m: list the unfinished tsv rows, ready for sbatch --array
# A crash that leaves a 0-byte file looks like empty. Needs gawk and bash 4.
check_all() {
	local list=$MC/scripts/src/datasets.tsv show_rows=0 G="" Y="" R="" O="" ds row tsv
	[[ $1 == -m ]] && { show_rows=1; shift; }
	local datasets=("$@")
	(( ${#datasets[@]} )) || mapfile -t datasets < <(grep -v '^#' "$list" | cut -f1)
	[[ -t 1 ]] && { G=$'\033[32m'; Y=$'\033[33m'; R=$'\033[31m'; O=$'\033[0m'; }

	# label;folders;depth under the folder;file name;removed from the file name to get the sample;what counts as a result
	local sm='_(RefSeq_20250528|gtdb-reps-rs220|gtdb-rs220-reps?)_gather[.]csv$'
	local mpa='_profile[.]txt$' kb='_bracken_S[.]MPA[.]TXT$' sy='_profile[.]tsv$'
	local rows=(
		"MPA22;MPA_db2022;2;*_profile.txt;$mpa;species"  "MPA23;MPA_db2023;2;*_profile.txt;$mpa;species"
		"MPA25;MPA_db2025;2;*_profile.txt;$mpa;species"
		"MOTU3;MOTUS;1;*_profile.txt;$mpa;motus"         "MOTU4;MOTUS4;1;*_profile.txt;$mpa;motus"
		"KB10;KB10;3;*_bracken_S.MPA.TXT;$kb;species"    "KB45;KB45;3;*_bracken_S.MPA.TXT;$kb;species"
		"KB90;KB90;3;*_bracken_S.MPA.TXT;$kb;species"    "KG10;KB10_GTDB;3;*_bracken_S.MPA.TXT;$kb;species"
		"KG45;KB45_GTDB;3;*_bracken_S.MPA.TXT;$kb;species" "KG90;KB90_GTDB;3;*_bracken_S.MPA.TXT;$kb;species"
		"SMref;SM_RefSeq_20250528;1;*_gather.csv;$sm;rows"
		"SM220;SM_gtdb-rs220-rep SM_gtdb-rs220-reps;1;*_gather.csv;$sm;rows"
		"SYref;SYLPH_RefSeq/per_sample;1;*_profile.tsv;$sy;rows" "SY220;SYLPH_gtdb_rs220/per_sample;1;*_profile.tsv;$sy;rows"
	)

	# One cell per row x dataset: "ok empty unfinished unfinished_rows"
	local -A cell total width
	local unfinished_list=()
	for ds in "${datasets[@]}"; do
		tsv=$(eval echo "$(awk -F'\t' -v ds="$ds" '$1 == ds { print $2 }' "$list")")
		if [[ ! -f $tsv ]]; then echo "No sample tsv for $ds (see scripts/src/datasets.tsv)"; return 1; fi
		total[$ds]=$(awk 'NF { n++ } END { print n + 0 }' "$tsv")
		width[$ds]=${#ds}
		for row in "${rows[@]}"; do
			local label folders depth name strip what dirs=() f
			IFS=';' read -r label folders depth name strip what <<< "$row"
			for f in $folders; do [[ -d $MC/data/$ds/$f ]] && dirs+=("$MC/data/$ds/$f"); done
			cell[$label,$ds]=$(
				{ (( ${#dirs[@]} )) && find "${dirs[@]}" -mindepth "$depth" -maxdepth "$depth" -type f -name "$name" -print0; } |
				xargs -0 -r awk -v strip="$strip" -v what="$what" '
					BEGINFILE { if (ERRNO) nextfile; id = FILENAME; sub(/.*\//, "", id); sub(strip, "", id); has = 0 }
					(what == "species" && /s__/ && !/^#/) ||
					(what == "rows" && FNR > 1) ||
					(what == "motus" && !/^#|^mOTU\t|unassigned/ && $NF + 0 > 0) { has = 1; nextfile }
					ENDFILE { print id "\t" has }' |
				awk -F'\t' 'NR == FNR { if (NF) id[++n] = $1; next }
					{ if ($2 || !($1 in seen)) seen[$1] = $2 }
					END {
						for (i = 1; i <= n; i++) {
							if (!(id[i] in seen)) { c++; rows = rows sep i; sep = "," }
							else if (seen[id[i]]) a++
							else b++
						}
						print a + 0, b + 0, c + 0, rows
					}' "$tsv" -
			)
			local ok empty left rows_left text
			read -r ok empty left rows_left <<< "${cell[$label,$ds]}"
			text="$ok/$empty/$left"
			(( ${#text} > ${width[$ds]} )) && width[$ds]=${#text}
			(( show_rows && left )) && unfinished_list+=("$ds $label: $rows_left")
		done
	done

	# Table, padded before colouring so that columns line up
	printf '%-6s' ""; for ds in "${datasets[@]}"; do printf '  %-*s' "${width[$ds]}" "$ds"; done; echo
	printf '%-6s' "n";  for ds in "${datasets[@]}"; do printf '  %-*s' "${width[$ds]}" "${total[$ds]}"; done; echo
	for row in "${rows[@]}"; do
		local label=${row%%;*} ok empty left rows_left colour
		printf '%-6s' "$label"
		for ds in "${datasets[@]}"; do
			read -r ok empty left rows_left <<< "${cell[$label,$ds]}"
			colour=$G; (( empty )) && colour=$Y; (( left )) && colour=$R
			printf '  %s%-*s%s' "$colour" "${width[$ds]}" "$ok/$empty/$left" "$O"
		done
		echo
	done
	echo "Cells: ok/empty/unfinished"
	(( show_rows )) && { (( ${#unfinished_list[@]} )) && printf '%s\n' "${unfinished_list[@]}" || echo "Nothing unfinished"; }
	return 0
}

# use_dataset <dataset>: set DATASET, DATASET_PATH, TSV and N_SAMPLES from scripts/src/datasets.tsv
use_dataset() {
	local list=$MC/scripts/src/datasets.tsv tsv
	if [[ -z $1 ]] || ! tsv=$(awk -F'\t' -v ds="$1" '$1 == ds { print $2; f = 1 } END { exit !f }' "$list"); then
		echo "Usage: use_dataset <dataset>, one of: $(grep -v '^#' "$list" | cut -f1 | xargs)"
		return 1
	fi
	export DATASET=$1
	export DATASET_PATH=$MC/data/$1
	export TSV=$(eval echo "$tsv")
	if [[ -f $TSV ]]; then
		export N_SAMPLES=$(wc -l < "$TSV" | tr -d ' ')
	else
		export N_SAMPLES=""
		echo "(sample tsv not there yet)"
	fi
	echo "DATASET=$DATASET  N_SAMPLES=$N_SAMPLES"
	echo "DATASET_PATH=$DATASET_PATH"
	echo "TSV=$TSV"
}

if [[ $(hostname) == "narval1" ]]; then

export ILAFORES=/scratch/ronj2303
export ILL_PIPELINES=/scratch/ronj2303/ILL_pipelines

export MC=$ILAFORES/MethodsComparison
export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
export MOSS=$MC/data/Moss
export PR19=$MC/data

fi
