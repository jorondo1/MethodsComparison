#!/bin/bash

#init
threads="12"
sample="false";
tsv="false";
kraken_db="/cvmfs/datahub.genap.ca/vhost34/def-ilafores/kraken2_dbs/k2_pluspfp_16gb_20210517"
bracken_readlen="150"
confidence="false"

SHORT_OPTS="ht:m:o:tsv:"
LONG_OPTS='help,kraken_db,bracken_readlen,confidence'

OPTS=$(getopt -o $SHORT_OPTS --long $LONG_OPTS -- "$@")

# make sure the params are entered correctly
if [ $? -ne 0 ];
then
    help_message;
    exit 1;
fi

# loop through input params
while true; do
    # echo $1
	case "$1" in
        -h | --help) help_message; exit 1; shift 1;;
        -t) threads=$2; shift 2;;
        -o) base_out=$2; shift 2;;
        -tsv) tsv=$2; shift 2;;
		--kraken_db) kraken_db=$2; shift 2;;
        --confidence) confidence=$2; shift 2;;
        --bracken_readlen) bracken_readlen=$2; shift 2;;
        --) help_message; exit 1; shift; break ;;
		*) break;;
	esac
done

singularity exec --writable-tmpfs -e \
	-B $ILL_PIPELINES:$ILL_PIPELINES \
	-B $ILAFORES:$ILAFORES \
	$ILL_PIPELINES/containers/kraken.2.1.2.sif bash -c "
while read -r sample fq1 fq2; do
outdir=${base_out}/${sample}
mkdir -p $outdir

kraken2 --memory-mapping \
	--confidence ${confidence} \
	--paired \
	--threads ${threads} \
	--db ${kraken_db} \
	--use-names \
	--output ${out_dir}/${sample}/${sample}_taxonomy_nt \
	--report ${out_dir}/${sample}/${sample}.kreport \
	${fq1} ${fq2}

rm ${out_dir}/${sample}/${sample}_taxonomy_nt

bracken \
    -d ${kraken_db} \
    -i ${out_dir}/${sample}/${sample}.kreport \
    -o ${out_dir}/${sample}/${sample}_bracken/${sample}_S.bracken \
    -w ${out_dir}/${sample}/${sample}_bracken/${sample}_bracken_S.kreport \
    -r $bracken_readlen 
    
done < "$tsv"
"