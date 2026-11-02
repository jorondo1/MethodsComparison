#!/usr/bin/env bash
# Fetch run-level ENA metadata for public metagenomic WGS runs (first public 2023+).
# One row per run; per-study summaries and ranking are done in rank_wgs_projects.R.
set -euo pipefail

root="$(cd "$(dirname "$0")/../../.." && pwd)"
out_dir="$root/Out/ena_candidates"
mkdir -p "$out_dir"

curl -sf -X POST https://www.ebi.ac.uk/ena/portal/api/search \
  --data-urlencode 'result=read_run' \
  --data-urlencode 'query=library_strategy="WGS" AND library_source="METAGENOMIC" AND tax_tree(408169) AND first_public>=2023-01-01 AND tax_id!=646099' \
  --data-urlencode 'fields=study_accession,sample_accession,scientific_name,sample_alias,sample_title,instrument_platform,library_layout,base_count,first_public,study_title' \
  --data-urlencode 'format=tsv' \
  --data-urlencode 'limit=0' \
  > "$out_dir/runs_2023.tsv"

echo "$(($(wc -l < "$out_dir/runs_2023.tsv") - 1)) runs written to $out_dir/runs_2023.tsv"
