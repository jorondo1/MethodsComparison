#!/bin/bash
# Pull from the HPC only the files 0_proc_data.R reads (+ the classification rate tables), into this repo.
# Usage, from the repo root:  bash scripts/sync_from_hpc.sh <user@host> [extra rsync options, e.g. --dry-run]
# Nothing is deleted locally, and local files newer than the HPC copy are kept (-u).
# Add a profiler here when 0_proc_data.R starts reading new files.
set -eo pipefail

HPC=${1:?usage: sync_from_hpc.sh <user@host> [rsync options]}; shift
REMOTE=/net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison

# First match wins; every parent directory of a wanted file must be included to be walked into.
rsync -avum "$@" \
	--include='/data/' \
	--include='/data/*/' \
	--include='/data/*/MOTUS*/' \
	--include='/data/*/MOTUS*/*_profile.txt' \
	--include='/data/*/SM_*/' \
	--include='/data/*/SM_*/*_gather.csv' \
	--include='/data/*/SM_*/*_lineages.csv' \
	--include='/data/*/SYLPH_*/' \
	--include='/data/*/SYLPH_*/merged_species.tsv' \
	--include='/data/*/SYLPH_*/classification_rate.tsv' \
	--include='/data/*/KB*/' \
	--include='/data/*/KB*/*/' \
	--include='/data/*/KB*/*/*_bracken/' \
	--include='/data/*/KB*/*/*_bracken/*_bracken_S.MPA.TXT' \
	--include='/data/*/MPA_*/' \
	--include='/data/*/MPA_*/*/' \
	--include='/data/*/MPA_*/*/*_profile.txt' \
	--include='/Out/' \
	--include='/Out/classification_rates/' \
	--include='/Out/classification_rates/*.tsv' \
	--exclude='*' \
	"$HPC:$REMOTE/" ./
