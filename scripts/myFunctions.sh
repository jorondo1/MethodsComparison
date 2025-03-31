#!/bin/bash

check_output() {
	if [ $# -lt 3 ]; then
	        echo "!!Usage: my_function <'database1 database2 ... databaseN'> <DATASETS> <filename_suffix>"
	        return 1  # Exit the function with a non-zero status
	fi
	filename_suffix="$3"
	DATASETS="$2"
	DATABASES="$1"

	for db in $DATABASES; do
		for ds in $DATASETS; do
			TSV_name="${ds}_TSV"
			TSV="${!TSV_name}"
			eval exp=\$NUM_$ds

			found=$(find $ds/*$db -type f -name "*$filename_suffix" -exec basename {} \; | sed "s/${filename_suffix}//" | sed "s/_${db}//")
			num=$(echo $found | wc -w)
			echo "$num $db output for $ds found, $exp expected."
			if [ "$num" -lt "$exp" ]; then 
			grep -nv -E "${found}" $TSV | cut -d: -f1| paste -s -d,
			fi
		done
	done
}

