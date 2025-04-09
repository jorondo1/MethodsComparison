#!/usr/bin/env bash

# Contact: Mike Lee (Mike.Lee@nasa.gov; github.com/AstrobioMike)

if [ "$#" != 0 ]; then
    printf "\n  Helper script to download all refseq complete genomes as of whatever today is.\n"
    printf "  See script for details. There are currently no guardrails or safety nets if a\n"
    printf "  download fails. So check the starting file count vs the total downloaded at end.\n\n"
    printf "  \tUsage:\n\t  bash refseq-complete-genome-dl.sh\n\n"

    printf "  \tContact:\n\t  Mike Lee (Mike.Lee@nasa.gov; github.com/AstrobioMike)\n\n"
    exit
fi

database="${1}"

# we can pull through https or ftp
# protocol="ftp"
protocol="https"

refseq_base_link="${protocol}://ftp.ncbi.nlm.nih.gov/refseq/release/${database}"
curr_date_marker=$(date +%d-%B-%Y)
refseq_html_file="refseq-${curr_date_marker}.html"
refseq_filenames_file="refseq-${curr_date_marker}-genome-files-${database}.txt"
genomes_dir="refseq-${curr_date_marker}-${database}-genomes"

mkdir -p ${genomes_dir}

# downloading html page (using this to get all the files we want to download)
curl -L -s -o ${refseq_html_file} ${refseq_base_link}

# parsing out genomic.fna.gz filenames (which are also their link suffixes)
grep "genomic.fna.gz" ${refseq_html_file} | cut -f 2 -d '"' > ${refseq_filenames_file}

# this is messy so that it works on darwin (mac) too
num_files=$(wc -l ${refseq_filenames_file} | sed 's/^ *//' | tr -s " " "\t" | cut -f 1)

printf "\n  We are beginning the download of ${num_files} files now...\n"
printf "  See you in a bit :)\n\n"

# downloading in parallel with xargs (num run in parallel is set with -P option)
xargs -I % -P 10 curl -L -s -O "${refseq_base_link}/%" < ${refseq_filenames_file}

mv *genomic.fna.gz ${genomes_dir}