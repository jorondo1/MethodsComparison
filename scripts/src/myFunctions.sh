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
