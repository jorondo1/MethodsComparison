#!/bin/bash
# Pull from the HPC only the files 0_proc_data.R reads (+ the classification rate tables), into this repo:
# the merged tables written by scripts/src/merge_tables.sh and sylph_merge.sh (run them on the HPC first), ~230 MB.
# Usage, from the repo root:  bash scripts/sync_from_hpc.sh <user@host> [extra rsync options, e.g. --dry-run]
# Nothing is deleted locally, and local files newer than the HPC copy are kept (-u).
# Add a profiler here when 0_proc_data.R starts reading new files.
# The other datasets' metadata (data/<ds>/raw/, P19 .rds) isn't synced: it's already local.
set -eo pipefail

HPC=${1:?usage: sync_from_hpc.sh <user@host> [rsync options]}; shift
REMOTE=/net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison

# First match wins; every parent directory of a wanted file must be included to be walked into.
rsync -avum "$@" \
	--include='/data/' \
	--include='/data/*/' \
	--include='/data/*/MOTUS*/' \
	--include='/data/*/MOTUS*/merged_profile.tsv' \
	--include='/data/*/KB*/' \
	--include='/data/*/KB*/merged_profile.tsv' \
	--include='/data/*/MPA_*/' \
	--include='/data/*/MPA_*/merged_profile.tsv' \
	--include='/data/*/SM_*/' \
	--include='/data/*/SM_*/merged_tax_summary.csv.gz' \
	--include='/data/*/SYLPH_*/' \
	--include='/data/*/SYLPH_*/merged_species.tsv' \
	--include='/data/*/SYLPH_*/classification_rate.tsv' \
	--include='/data/Maple/raw/' \
	--include='/data/Maple/raw/metadata.tsv' \
	--include='/Out/' \
	--include='/Out/classification_rates/' \
	--include='/Out/classification_rates/*.tsv' \
	--exclude='*' \
	"$HPC:$REMOTE/" ./
