#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH -D /net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison
#SBATCH -o /net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison/logs/sylph_merge-%j.slurm.out
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A def-ilafores
#SBATCH -J sylph_merge

# Merge the per-sample outputs of sylph_SLURM.sh into tables. Run once all sylph_SLURM.sh array tasks are done
# (from $MC, or sbatch --dependency=afterany:<array job id>).
# Usage: sylph_merge.sh <dataset> <db name> <sample tsv>
# Writes to data/<dataset>/SYLPH_<db name>/:
#   counts.tsv           : estimated read counts, genomes (rows) x samples (columns), corrected (see count correction below)
#   taxonomy.tsv         : genome -> Domain, Phylum, Class, Order, Family, Genus, Species (same rows as counts.tsv)
#   merged_all_ranks.tsv : sylph-tax merge of the per-sample profiles, corrected read counts at every rank
#   classification_rate.tsv : path, % classified, % unclassified, rate (no header, like the other tools' rate files)
set -eo pipefail

export DB_NAME="${2}"
export SAM_LIST="${3}"
export OUT_DIR=${PWD}/data/"${1}"/SYLPH_${DB_NAME}
export SAMPLE_DIR=${OUT_DIR}/per_sample

# sylph + sylph-tax: conda env under $ILAFORES, so every node sees it (unlike ~/miniconda3, missing on the old nodes).
# Putting its bin/ first in PATH is enough for these two command-line tools: no conda needed in the job.
export SYLPH_ENV=${SYLPH_ENV:-$ILAFORES/programs/conda_envs/sylph}
export PATH=$SYLPH_ENV/bin:$PATH
command -v sylph-tax > /dev/null || { echo "sylph-tax not found on $(hostname) (SYLPH_ENV=$SYLPH_ENV)"; exit 1; }

# Check every sample of the tsv: no profile = sylph_SLURM.sh didn't finish (rerun it), no .sylphmpa = no genome detected
missing=0
for s in $(cut -f1 $SAM_LIST); do
	if [[ ! -f $SAMPLE_DIR/${s}_profile.tsv ]]; then
		echo "$s: no profile, rerun sylph_SLURM.sh for it"; missing=$((missing + 1))
	elif [[ ! -f $SAMPLE_DIR/${s}.sylphmpa ]]; then
		echo "$s: no genome detected, left out of the tables"
	fi
done
if [[ $missing -gt 0 ]]; then
	echo "$missing samples unfinished, not merging"
	exit 1
fi

# Classification rate = % of reads detected in database, as logged by sylph profile (-u, forced by --estimate-read-counts):
# bases of the sample covered by the detected genomes / all bases of the sample
echo "Classification rates..."
rm -f $OUT_DIR/classification_rate.tsv
no_rate=0
for s in $(cut -f1 $SAM_LIST); do
	log=$SAMPLE_DIR/${s}_profile.log
	path=./data/${1}/SYLPH_${DB_NAME}/per_sample/${s}_profile.log
	pct=$( [[ -f $log ]] && grep -o '[0-9.]*% of reads detected in database' $log | tail -1 | cut -d% -f1 || true )
	if [[ -n $pct ]]; then
		awk -v path=$path -v pct=$pct 'BEGIN { printf "%s\t%.2f\t%.2f\t%.5f\n", path, pct, 100 - pct, pct / 100 }' >> $OUT_DIR/classification_rate.tsv
	else
		printf "%s\tNA\tNA\tNA\n" $path >> $OUT_DIR/classification_rate.tsv
		# the rate is needed to correct the counts of every sample that made it into the tables
		if [[ -f $SAMPLE_DIR/${s}.sylphmpa ]]; then echo "$s: no classification rate in $log"; no_rate=$((no_rate + 1)); fi
	fi
done
if [[ $no_rate -gt 0 ]]; then
	echo "$no_rate samples without classification rate, can't correct their counts: not merging"
	echo "(profiled before _profile.log was kept? recover the logs from the slurm logs, see MethodsComparison_Pipeline.sh)"
	exit 1
fi

echo "Merging $(ls $SAMPLE_DIR/*.sylphmpa | wc -l) samples into $OUT_DIR..."
sylph-tax merge $SAMPLE_DIR/*.sylphmpa --column sequence_abundance -o $OUT_DIR/merged_all_ranks.raw.tsv

# Count correction: with --estimate-read-counts, sylph computes each genome's count as
# coverage x genome size / read length x rate, and coverage x genome size / read length is already the read count.
# The extra x rate makes counts sum to rate^2 x reads instead of rate x reads (checked on Maple 1_1_WGS_20250513:
# counts/reads 0.005, rate 0.076, rate^2 0.006), so each sample's counts are divided by its rate.
awk -F'\t' -v OFS='\t' '
	NR == FNR { s = $1; sub(/.*\//, "", s); sub(/_profile\.log$/, "", s); rate[s] = $4; next }
	FNR == 1 { for (i = 2; i <= NF; i++) r[i] = rate[$i]; print; next }
	{ for (i = 2; i <= NF; i++) $i = (r[i] > 0) ? sprintf("%.0f", $i / r[i]) : 0; print }
	' $OUT_DIR/classification_rate.tsv $OUT_DIR/merged_all_ranks.raw.tsv > $OUT_DIR/merged_all_ranks.tsv
rm $OUT_DIR/merged_all_ranks.raw.tsv

# genome rows (t__) only: clade_name = d__X|p__X|...|s__X|t__<genome>, or NO_TAXONOMY|t__<genome>
awk -F'\t' -v counts=$OUT_DIR/counts.tsv -v tax=$OUT_DIR/taxonomy.tsv '
	BEGIN { split("d p c o f g s", prefix, " "); split("Domain Phylum Class Order Family Genus Species", rank, " ") }
	NR == 1 {
		printf "genome" > counts; for (i = 2; i <= NF; i++) printf "\t%s", $i > counts; printf "\n" > counts
		printf "genome" > tax; for (r = 1; r <= 7; r++) printf "\t%s", rank[r] > tax; printf "\n" > tax
		next
	}
	$1 ~ /(^|\|)t__/ {
		n = split($1, lev, "|"); genome = substr(lev[n], 4)
		printf "%s", genome > counts; for (i = 2; i <= NF; i++) printf "\t%s", $i > counts; printf "\n" > counts
		delete name
		for (j = 1; j < n; j++) name[substr(lev[j], 1, 1)] = substr(lev[j], 4)
		printf "%s", genome > tax; for (r = 1; r <= 7; r++) printf "\t%s", name[prefix[r]] > tax; printf "\n" > tax
	}' $OUT_DIR/merged_all_ranks.tsv

echo "$(head -1 $OUT_DIR/counts.tsv | awk -F'\t' '{print NF-1}') samples, $(tail -n+2 $OUT_DIR/counts.tsv | wc -l) genomes"
echo "Merged tables written to $OUT_DIR"
