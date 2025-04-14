#!/bin/bash

# Verify output for all datasets
check_output() {
	if [ $# -lt 3 ]; then
	        echo "!!Usage: my_function <'database1 database2 ... databaseN'> <filename_suffix>"
	        echo "Depends on variables exported by DATASET function."
	        return 1  # Exit the function with a non-zero status
	fi
	
    if [[ -z "${!DATASET+x}" || -z "${!TSV+x}" ]]; then
		echo "\$DATASET not set. Run dataset <DATASET> <TSV> to declare variables."
		return 1
	fi
	
	DATABASES="$1"
	filename_suffix="$2"

	for db in $DATABASES; do
		
		found=$(find $DATASET/*$db -type f -name "*$filename_suffix" -exec basename {} \; | sed "s/${filename_suffix}//" | sed "s/_${db}//")
		
		num=$(echo $found | wc -w)
		echo "$num $db output for $DATASET found, $N_SAMPLES expected."
		
		if [ "$num" -lt "$exp" ]; then 
			grep -nv -E "${found}" $TSV | cut -d: -f1| paste -s -d,
		fi
	done
}

# Function to export variable names
dataset_variables() {
    declare -g "DATASET"="$1"
    declare -g "TSV"="$2"
    declare -g "N_SAMPLES"=$(wc -l < "$2")
    echo "\$DATASET evaluates to $1"
    echo "\$TSV evaluates to $TSV"
    echo "\$N_SAMPLES evaluates to $N_SAMPLES"
}


