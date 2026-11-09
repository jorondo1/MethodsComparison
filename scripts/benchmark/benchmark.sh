# Resource benchmark: runtime and memory per tool x database x dataset, on one node type (c-3t).
# Copy-paste collection: running or sourcing this file does nothing. Everything below goes through sbatch.
# Output: Out/benchmark/ (never data/). Figure: scripts/10_benchmark_plot.R
return 0 2>/dev/null; [[ $0 == "${BASH_SOURCE[0]}" ]] && { echo "Copy-paste collection: paste the commands, don't run this file"; exit 1; }

export MC=$ILAFORES/analysis/MethodsComparison
cd $MC

B=$MC/Out/benchmark # output
BM=$MC/scripts/benchmark/src
SRC=$MC/scripts/src
DB_MPA=$FAST/metaphlan4_db
DB_SM=$ILAFORES/ref_dbs/sourmash_db
DB_SYLPH=$ILAFORES/ref_dbs/sylph_db
DB_K2=$ILAFORES/ref_dbs/kraken2_dbs

N=5 # samples per dataset
MIN_MB=10 # skip near-empty samples (paired_1 size): they say nothing about runtime, and some fail

################################################################
# Pick samples ########
################################################################

# samples.tsv: dataset, sample, 4 read files, paired_1 size. N rows per dataset, evenly spaced by read file size,
# leaving out samples with a missing read file or under MIN_MB.
# tools.tsv: the same rows, as the 5 columns the tool scripts read. Array index = row of these files.
mkdir -p $B/logs
grep -v '^#' scripts/src/datasets.tsv | while IFS=$'\t' read -r ds tsv; do
	tsv=$(eval echo "$tsv")
	while IFS=$'\t' read -r sample p1 p2 u1 u2 _; do
		if [[ ! ( -s $p1 && -s $p2 && -s $u1 && -s $u2 ) ]]; then
			echo "Skipped $ds $sample: missing read file" >&2
		elif (( $(stat -c %s "$p1") >= MIN_MB * 1024 * 1024 )); then
			printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(stat -c %s "$p1")" "$sample" "$p1" "$p2" "$u1" "$u2"
		fi
	done < "$tsv" | sort -n | awk -F'\t' -v OFS='\t' -v n=$N -v ds=$ds '
		{ row[NR] = $0 }
		END {
			for (i = 1; i <= n && i <= NR; i++) {
				k = (NR <= n) ? i : int(1 + (i - 1) * (NR - 1) / (n - 1) + 0.5)
				split(row[k], f, "\t")
				print ds, f[2], f[3], f[4], f[5], f[6], f[1]
			}
		}'
done > $B/samples.tsv
cut -f2-6 $B/samples.tsv > $B/tools.tsv

A=1-$(wc -l < $B/tools.tsv) # array rows. Pilot: A=1
K=$(wc -l < $B/tools.tsv)   # Kraken rows.  Pilot: K=1

# submit <label> <sbatch args>: runs on c-3t, logs to $B/logs, adds label + job id to jobs.tsv, sets $LAST
submit() {
	local label=$1; shift
	LAST=$(sbatch --parsable -p c-3t -D $MC -J bm_$label --mail-type=NONE \
		-o $B/logs/$label-%A_%a.out -e $B/logs/$label-%A_%a.err "$@") || return 1
	LAST=${LAST%%;*}
	printf '%s\t%s\n' "$label" "$LAST" >> $B/jobs.tsv
	echo "$label $LAST"
}

################################################################
# Bases per sample ########
################################################################

submit bases $BM/benchmark_bases.sh $B # writes $B/bases.tsv

################################################################
# Sketches, as their own rows ########
################################################################

submit sylph_sketch    --array=$A%10 -n 2 --mem=16G $BM/benchmark_sketch.sh sylph $B/tools.tsv $B
SK_SYLPH=$LAST
submit sourmash_sketch --array=$A%10 -n 1 --mem=16G $BM/benchmark_sketch.sh sourmash $B/tools.tsv $B
SK_SM=$LAST

################################################################
# MetaPhlAn ########
################################################################

# ILL_pipelines must be on kneaddata_gzip, as in pipeline_2_profiling.sh
metaphlan="bash $ILL_PIPELINES/generateslurm_taxonomic_abundance.metaphlan.4.2.4.sh \
	--slurm_log $B/logs --slurm_walltime 48:00:00 --slurm_threads 12 --slurm_mem 100G --sample_tsv $B/tools.tsv"
$metaphlan --db $DB_MPA/2022/mpa_vOct22_CHOCOPhlAnSGB_202212 --out $B/MPA_db2022
$metaphlan --db $DB_MPA/2023/mpa_vJun23_CHOCOPhlAnSGB_202307 --out $B/MPA_db2023
$metaphlan --db $DB_MPA/2025/mpa_vJan25_CHOCOPhlAnSGB_202503 --out $B/MPA_db2025

submit MPA22 --array=$A%2 $B/MPA_db2022/metaphlan.slurm.sh
submit MPA23 --array=$A%2 $B/MPA_db2023/metaphlan.slurm.sh
submit MPA25 --array=$A%2 $B/MPA_db2025/metaphlan.slurm.sh

################################################################
# mOTUs ########
################################################################

submit MOTU3 --array=$A%3 -n 6  --mem=60G  $SRC/motus3_SLURM.sh $B $B/tools.tsv
submit MOTU4 --array=$A%2 -n 12 --mem=100G $SRC/motus4_SLURM.sh $B $B/tools.tsv

################################################################
# sylph, sourmash ########
################################################################

# afterany: a failed sketch then costs that sample only, not the whole array
# ../Out/benchmark: these scripts write to data/<dataset>/, so this lands in $B, next to the sketches above
DS=../Out/benchmark

submit SYref --array=$A%5 -n 2 --mem=64G --dependency=afterany:$SK_SYLPH $SRC/sylph_SLURM.sh \
	$DS $B/tools.tsv $DB_SYLPH/RefSeq_20250528_c200.syldb $DB_SYLPH/RefSeq_20250528_taxonomy.tsv.gz RefSeq
submit SY220 --array=$A%5 -n 2 --mem=64G --dependency=afterany:$SK_SYLPH $SRC/sylph_SLURM.sh \
	$DS $B/tools.tsv $DB_SYLPH/gtdb-r220-c200-dbv1.syldb $DB_SYLPH/gtdb_r220_metadata.tsv.gz gtdb_rs220

submit SMref --array=$A%5 -n 1  --mem=64G --dependency=afterany:$SK_SM $SRC/gather_SLURM_fast.sh \
	$DS $B/tools.tsv $DB_SM/RefSeq_20250528.k31.sbt.zip
submit SM220 --array=$A%2 -n 24 --mem=64G --dependency=afterany:$SK_SM $SRC/gather_SLURM_fast.sh \
	$DS $B/tools.tsv $DB_SM/gtdb-rs220-reps-k31.zip

################################################################
# Kraken + Bracken ########
################################################################

# One job per index. Copies it to /dev/shm, then one srun step per sample x confidence, one after the other.
# args: index, step prefix (KB or KG), output folder suffix, tsv, number of rows, output
submit KB --time=48:00:00 -N 1 -c 24 --mem=200G $BM/benchmark_kraken.sh $DB_K2/k2_standard_20241228 KB "" $B/tools.tsv $K $B
submit KG --time=48:00:00 -N 1 -c 24 --mem=700G $BM/benchmark_kraken.sh $DB_K2/k2_gtdb_genome_reps_20241109 KG _GTDB $B/tools.tsv $K $B

################################################################
# Once everything has finished ########
################################################################

check_jobs() { squeue -u $USER -h -o '%j %T' | grep '^bm_' | sort | uniq -c; } # empty = finished
check_jobs

bash $BM/benchmark_gather.sh $B # writes $B/sacct.tsv

module load StdEnv/2023 r/4.4.0
VROOM_THREADS=2 Rscript scripts/10_benchmark_plot.R $B
