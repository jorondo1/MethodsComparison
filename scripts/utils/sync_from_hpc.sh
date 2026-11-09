#!/bin/bash
# Pull from the HPC the inputs of the R scripts, into this repo:
#   - the merged tables written by scripts/src/merge_tables.sh and scripts/utils/sylph_merge.sh (run them on the HPC first), ~230 MB;
#   - the sample data read by scripts/utils/sample_data.R: the top-level .tsv/.csv/.txt of data/<ds>/raw/ (~1 MB in all;
#     the read folders inside raw/ are not walked into) and data/<ds>/*.rds (P19 phyloseq objects);
#   - each dataset's sample tsv (preproc/preprocessed_reads.sample.tsv, or raw/ for Moss), read by the 0_run_all.R checks;
#   - the classification rate tables (Out/classification_rates/*.tsv, every tool);
#   - the phyloseq lists and diversity results scripts/0_run_all.R assembles on the HPC (Out/_Rdata/*.RDS, listed by
#     name so the stale 2025 leftovers there are not pulled; the per-dataset cache, Out/_Rdata/cache/, stays on the HPC).
# Usage, from the repo root:  bash scripts/utils/sync_from_hpc.sh <user@host> [extra rsync options, e.g. --dry-run]
# Nothing is deleted locally, and local files newer than the HPC copy are kept (-u).
# When an R script starts reading a new kind of file, or 0_run_all.R writes a new RDS, add it here.
set -eo pipefail

HPC=${1:?usage: sync_from_hpc.sh <user@host> [rsync options]}; shift
REMOTE=/net/nfs-ip34/jbod2/def-ilafores/analysis/MethodsComparison

# First match wins; every parent directory of a wanted file must be included to be walked into.
rsync -avum "$@" \
	--exclude='._*' \
	--include='/data/' \
	--include='/data/*/' \
	--include='/data/*/*.rds' \
	--include='/data/*/raw/' \
	--include='/data/*/raw/*.tsv' \
	--include='/data/*/raw/*.csv' \
	--include='/data/*/raw/*.txt' \
	--include='/data/*/preproc/' \
	--include='/data/*/preproc/preprocessed_reads.sample.tsv' \
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
	--include='/data/*/SYLPH_*/merged_species_tax.tsv' \
	--include='/data/*/SYLPH_*/classification_rate.tsv' \
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
