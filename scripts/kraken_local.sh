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

#init
threads="12"
tsv="false";
kraken_db="$ILAFORES/ref_dbs/kraken2_dbs/k2_standard_20241228"
bracken_readlen="150"
confidence="0.05"

SHORT_OPTS="h:t:o:tsv:"
LONG_OPTS='help,kraken_db:,bracken_readlen:,confidence:'

OPTS=$(getopt -o $SHORT_OPTS --long $LONG_OPTS -- "$@")

# make sure the params are entered correctly
if [ $? -ne 0 ];
then
    help_message;
    exit 1;
fi

# loop through input params
eval set -- "$OPTS"
while true; do
    case "$1" in
        -h | --help) help_message; exit 0;;
        -t) threads=$2; shift 2;;
        -o) out_dir=$2; shift 2;;
        -tsv) tsv=$2; shift 2;;
        --kraken_db) kraken_db=$2; shift 2;;
        --confidence) confidence=$2; shift 2;;
        --bracken_readlen) bracken_readlen=$2; shift 2;;
        --) shift; break;;
        *) echo "Invalid option: $1"; help_message; exit 1;;
    esac
done

# Check required parameters
if [ "$tsv" = "false" ]; then
    echo "Error: TSV file not provided"
    help_message
    exit 1
fi

if [ ! -f "$tsv" ]; then
    echo "Error: TSV file $tsv does not exist"
    exit 1
fi

singularity exec --writable-tmpfs -e \
    -B $ILL_PIPELINES:$ILL_PIPELINES \
    -B $ILAFORES:$ILAFORES \
    $ILL_PIPELINES/containers/kraken.2.1.2.sif bash -c "
while read -r sample fq1 fq2; do
    out_dir=\"${out_dir}/\${sample}\"
    mkdir -p \"\$out_dir\"

    kraken2 --memory-mapping \
        --confidence ${confidence} \
        --paired \
        --threads ${threads} \
        --db \"${kraken_db}\" \
        --use-names \
        --output \"\${out_dir}/\${sample}_taxonomy_nt\" \
        --report \"\${out_dir}/\${sample}.kreport\" \
        \"\${fq1}\" \"\${fq2}\"

    rm \"\${out_dir}/\${sample}_taxonomy_nt\"

    mkdir -p \"\${out_dir}/\${sample}_bracken\"
    bracken \
        -d \"${kraken_db}\" \
        -i \"\${out_dir}/\${sample}.kreport\" \
        -o \"\${out_dir}/\${sample}_bracken/\${sample}_S.bracken\" \
        -w \"\${out_dir}/\${sample}_bracken/\${sample}_bracken_S.kreport\" \
        -r $bracken_readlen 
    
done < \"$tsv\"
"