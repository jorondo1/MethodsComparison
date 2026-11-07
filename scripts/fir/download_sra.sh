#!/bin/bash
# Download SRA runs from NCBI's public S3 copy (much faster than ENA's FTP). Needs internet: run on a login node, in tmux.
# Usage: bash download_sra.sh <filereport.tsv> [parallel=4]
#   filereport: ENA filereport with a header and run_accession in column 1. Files go to ./sra/<run>
# A file is done when its size equals the server's Content-Length. Partial files resume; finished runs are skipped.
# Measured on Fir login2: ~40 MB/s per stream, linear up to 4 streams.
set -uo pipefail

FR=${1:?filereport.tsv}
P=${2:-4}
mkdir -p sra

get() {
	local run=$1 url=https://sra-pub-run-odp.s3.amazonaws.com/sra/$1/$1 out=sra/$1
	local len size try
	len=$(curl -sfI --retry 3 "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
	[[ -n $len ]] || { printf '%s\tno_url\n' "$run"; return; }
	for try in 1 2 3; do
		size=$(stat -c %s "$out" 2>/dev/null || echo 0)
		(( size == len )) && { printf '%s\tok\t%s\n' "$run" "$len"; return; }
		(( size > len )) && rm -f "$out"
		curl -sf --retry 3 -C - -o "$out" "$url"
	done
	size=$(stat -c %s "$out" 2>/dev/null || echo 0)
	(( size == len )) && printf '%s\tok\t%s\n' "$run" "$len" || printf '%s\tincomplete\t%s/%s\n' "$run" "$size" "$len"
}
export -f get

tail -n +2 "$FR" | cut -f1 | xargs -P "$P" -I{} bash -c 'get {}' | tee download_log.tsv
echo "Done: $(grep -c $'\tok\t' download_log.tsv) ok of $(( $(wc -l < "$FR") - 1 )). Rerun to retry the rest."
