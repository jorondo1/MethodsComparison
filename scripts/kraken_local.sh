#!/bin/bash

help_message () {
    echo ""
    echo "Options:"

    echo ""
    echo "	  -tsv STR  path to tsv with 3 columns: sample_name path_to_fastq1 path_to_fastq2"
    echo "    -o STR    path to output dir"
    echo "    -t    # of threads (default 12)"
    echo "    --kraken_db    kraken2 database path (default /cvmfs/datahub.genap.ca/vhost34/def-ilafores/kraken2_dbs/k2_pluspfp_16gb_20210517)"
    echo "    --bracken_readlen    bracken read length option (default 150)"
    echo "    --confidence    kraken confidence level to reduce false-positive rate (default 0.05)"
    echo "    -tsv STR    input TSV file with sample names and paths"

    echo ""
    echo "  -h --help    Display help"

    echo "";
}

# Initialize with empty values (no defaults for required params)
threads=""
tsv=""
kraken_db=""
bracken_readlen="150"
confidence="0.05"
out_dir=""

SHORT_OPTS="h:t:o:"
LONG_OPTS='help:,threads:,output:,kraken_db:,bracken_readlen:,confidence:,tsv:'

OPTS=$(getopt -o $SHORT_OPTS --long $LONG_OPTS -- "$@") || { help_message; exit 1; }

eval set -- "$OPTS"
while true; do
    case "$1" in
        -h | --help) help_message; exit 0;;
        -t | --threads) threads="$2"; shift 2;;
        -o | --output) out_dir="$2"; shift 2;;
        --tsv) tsv="$2"; shift 2;;
        --kraken_db) kraken_db="$2"; shift 2;;
        --confidence) confidence="$2"; shift 2;;
        --bracken_readlen) bracken_readlen="$2"; shift 2;;
        --) shift; break;;
        *) echo "Invalid option: $1"; help_message; exit 1;;
    esac
done

# Verify all required parameters
missing=()
[[ -z "$tsv" ]] && missing+=("--tsv")
[[ -z "$out_dir" ]] && missing+=("--output")
[[ -z "$kraken_db" ]] && missing+=("--kraken_db")

if (( ${#missing[@]} > 0 )); then
    echo "ERROR: Missing required options: ${missing[*]}"
    help_message
    exit 1
fi

# Create output directory if needed
mkdir -p "$out_dir"

# Start logging
exec > >(tee -a "${out_dir}/kraken_wrapper.log") 2>&1

# Copying database to shared memory
#rsync -avr $kraken_db /dev/shm

while IFS=$'\t' read -r sample fq1 fq2 _; do
    total_start=$(date +%s)
    iter_start=$(date +%s)
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Starting sample \${sample}..."
    
    out_dir="${out_dir}/${sample}"
    mkdir -p "$out_dir"

    if [[ -f "${out_dir}/${sample}_bracken/${sample}_bracken_S.MPA.TXT" ]]; then
        echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Skipping ${sample} - outputs already exist" 
        continue
    fi

singularity exec --writable-tmpfs -e \
    -B $ILL_PIPELINES:$ILL_PIPELINES \
    -B /dev:/dev \
    -B $ILAFORES:$ILAFORES \
    $ILAFORES/programs/ILL_pipelines/containers/kraken.2.1.2.sif bash -c "
    echo \"Running with ${threads} threads\"
    # Kraken2 with memory mapping
    
    last_core=\$(( ${threads} - 1 ))
    
    # Execute Kraken2 with taskset
    taskset -c 0-\${last_core} kraken2 --memory-mapping \\
        --confidence ${confidence} \\
        --paired \\
        --threads \"${threads}\" \\
        --db \"${kraken_db}\" \\
        --use-names \\
        --output \"${out_dir}/${sample}_taxonomy_nt\" \\
        --report \"${out_dir}/${sample}.kreport\" \\
        \"${fq1}\" \"${fq2}\"

    # Bracken processing
    mkdir -p \"${out_dir}/${sample}_bracken\"
    bracken \\
        -d \"${kraken_db}\" \\
        -i \"${out_dir}/${sample}.kreport\" \\
        -o \"${out_dir}/${sample}_bracken/${sample}_bracken_S.MPA.TXT\" \\
        -w \"${out_dir}/${sample}_bracken/${sample}_bracken_S.kreport\" \\
        -r $bracken_readlen

    #rm \"${out_dir}/${sample}_taxonomy_nt\"
"
    iter_end=$(date +%s)
    iter_time=$((iter_end - iter_start))
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Completed ${sample} in ${iter_time} seconds"
    
done < "$tsv"

total_end=$(date +%s)
total_time=$((total_end - total_start))
echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] All samples processed in ${total_time} seconds"
