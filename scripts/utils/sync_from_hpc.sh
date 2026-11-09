#!/bin/bash
# Pull from the HPC the inputs of the figure scripts, into this repo. The pipeline (merge_tables.sh, sylph_merge.sh,
# scripts/0_run_all.R) runs on the HPC only, so its inputs (merged tables, sample data, sample tsvs) are not synced:
#   - the classification rate tables (Out/classification_rates/*.tsv, every tool);
#   - the phyloseq lists and diversity results scripts/0_run_all.R assembles (Out/_Rdata/*.RDS, listed by name so the
#     stale 2025 leftovers there are not pulled; the per-dataset cache, Out/_Rdata/cache/, stays on the HPC).
# scripts/src/datasets.tsv, also read by 8_Identifiability.R, comes with git.
# Usage, from the repo root:  bash scripts/utils/sync_from_hpc.sh <user@host> [extra rsync options, e.g. --dry-run]
# Nothing is deleted locally, and local files newer than the HPC copy are kept (-u).
# When 0_run_all.R writes a new RDS, or a figure script reads a new kind of file, add it here.
set -eo pipefail

HPC=${1:?usage: sync_from_hpc.sh <user@host> [rsync options]}; shift
REMOTE=/net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison

# First match wins; every parent directory of a wanted file must be included to be walked into.
rsync -avum "$@" \
	--include='/Out/' \
	--include='/Out/classification_rates/' \
	--include='/Out/classification_rates/*.tsv' \
	--include='/Out/_Rdata/' \
	--include='/Out/_Rdata/ps_raw.ls.RDS' \
	--include='/Out/_Rdata/ps_rare.ls.RDS' \
	--include='/Out/_Rdata/alpha_div.RDS' \
	--include='/Out/_Rdata/pcoa_noVST.ls.RDS' \
	--include='/Out/_Rdata/pairwise_dist.RDS' \
	--include='/Out/_Rdata/procrustes.RDS' \
	--include='/Out/_Rdata/permanova.ds.RDS' \
	--include='/Out/_Rdata/dataset_names_n_lookup_alphadiv.RDS' \
	--exclude='*' \
	"$HPC:$REMOTE/" ./
