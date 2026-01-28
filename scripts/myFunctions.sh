#!/bin/bash

# Verify output for all datasets
check_output() {
	if [ $# -lt 2 ]; then
	        echo "!!Usage: my_function <'database1 database2 ... databaseN'> <filename_suffix>"
	        echo "Depends on variables exported by DATASET function."
	        return 1  # Exit the function with a non-zero status
	fi
	
    if [[ -z "$DATASET" || -z "$TSV" ]]; then
		echo "\$DATASET or \$TSV not set. Run dataset_variables <DATASET> <TSV> to declare variables."
		return 1
	fi
	
	DATABASES="$1"
	DATASET="$2"
	filename_suffix="$3"

	for db in $DATABASES; do
		
		found=$(find $DATASET/*$db -type f -name "*$filename_suffix" -exec basename {} \; | sed "s/${filename_suffix}//" | sed "s/_${db}//")
		
		num_found=$(echo $found | wc -w)
		echo "$num_found $db output for $DATASET found, $N_SAMPLES expected."
		
		if [ "$num_found" -lt "$N_SAMPLES" ]; then 
			export FOUND=$(grep -nv -E "${found}" "$TSV" | cut -d: -f1| paste -s -d,)
			echo $FOUND
		fi
	done
}

# Function to export variable names
dataset_variables() {
    declare -g "DATASET"="$1"
    declare -g "TSV"="$2"
    declare -g "N_SAMPLES"=$(wc -l < "$2")
    declare -g "DATASET_PATH"="$MC/data/$DATASET"
    echo "\$DATASET evaluates to $1"
    echo "\$TSV evaluates to $TSV"
    echo "\$N_SAMPLES evaluates to $N_SAMPLES"
    echo "\$DATASET_PATH evaluates to $DATASET_PATH"
}

if [[ $(hostname) == "narval1" ]]; then

export ILAFORES=/scratch/ronj2303
export ILL_PIPELINES=/scratch/ronj2303/ILL_pipelines

export MC=$ILAFORES/MethodsComparison
export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
export MOSS=$MC/data/Moss
export PR19=$MC/data

fi
