#!/bin/bash
# Convert one downloaded .sra (sra/<run>) to <run>_1/_2.fastq.gz in the current directory. No internet needed.
# Slurm array, one run per task (row N of the filereport, header excluded). Submit from the folder that holds sra/:
#   sbatch -A def-ilafores -c 4 --mem 8G --time 1:00:00 --array=1-N -o logs/convert-%A_%a.out convert_sra.sh filereport.tsv
# Checks the _1 read count against the filereport's read_count, logs to download_status.tsv, deletes the .sra when ok.
set -euo pipefail
module load StdEnv/2023 sra-toolkit/3.0.9

FR=$(realpath "${1:?filereport.tsv}")
I=${SLURM_ARRAY_TASK_ID:?run as an array}
T=${SLURM_CPUS_PER_TASK:-4}
TMP=${SLURM_TMPDIR:-/tmp}/fqd_$I

col=$(head -1 "$FR" | tr '\t' '\n' | grep -nx read_count | cut -d: -f1)
read -r run expected < <(awk -F'\t' -v n=$((I + 1)) -v c="$col" 'NR==n{print $1, $c}' "$FR")
[[ -n $run ]] || { echo "No row $I in $FR"; exit 1; }
echo "$(date) task $I: $run, expected $expected reads"

if [[ -s ${run}_1.fastq.gz && -s ${run}_2.fastq.gz ]] && grep -q "^${run}"$'\tok' download_status.tsv 2>/dev/null; then
	echo "Already converted"; exit 0
fi
[[ -s sra/$run/$run || -s sra/$run ]] || { echo "Missing sra/$run"; exit 1; }
SRA=sra/$run; [[ -d $SRA ]] && SRA=sra/$run/$run

mkdir -p "$TMP"
fasterq-dump --split-3 -e "$T" -t "$TMP" -O "$TMP" \
	--seq-defline '@$ac.$si/$ri' --qual-defline '+' "$(realpath "$SRA")"

for r in 1 2; do
	pigz -p "$T" -c "$TMP/${run}_$r.fastq" > "${run}_$r.fastq.gz.tmp"
done
n=$(( $(wc -l < "$TMP/${run}_1.fastq") / 4 ))
rm -rf "$TMP"

if [[ $n == "$expected" ]]; then
	for r in 1 2; do mv "${run}_$r.fastq.gz.tmp" "${run}_$r.fastq.gz"; done
	printf '%s\tok\t%s\t%s\n' "$run" "$n" "$expected" >> download_status.tsv
	rm -f "sra/$run"; rm -rf "sra/$run"
	echo "$(date) ok: $n reads"
else
	printf '%s\tbad\t%s\t%s\n' "$run" "$n" "$expected" >> download_status.tsv
	echo "$(date) bad: $n reads, expected $expected. Kept sra/$run"; exit 1
fi
