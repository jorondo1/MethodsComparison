#!/usr/bin/env python3
"""Multi-sample `sourmash tax metagenome -F csv_summary`.

The sourmash CLI writes csv_summary for one gather CSV at a time ("Can only output
multi-query result formats"), and reloads the lineage CSVs at each call. This script runs
the same sourmash functions as `sourmash tax metagenome` (4.8.11, sourmash/tax/__main__.py)
without that check: lineages are loaded once and tax_utils.write_summary() writes every
query into one csv_summary file. The output equals the per-sample csv_summary files
concatenated (one header).

Run it with the python of the sourmash container:
  singularity exec -B $ILAFORES sourmash.4.8.11.sif python3 sourmash_tax_summary.py \
      -t lineages.csv [lineages2.csv ...] -o out.csv [--keep-full-identifiers] \
      [--keep-identifier-versions] gather1.csv gather2.csv ...
  (or --from-file list.txt instead of the gather CSVs)
"""
import argparse
import sys

from sourmash.tax import tax_utils
from sourmash.tax.tax_utils import MultiLineageDB
from sourmash.sourmash_args import FileOutputCSV


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("gather_csv", nargs="*")
    p.add_argument("--from-file", help="text file listing gather CSVs, one per line")
    p.add_argument("-t", "--taxonomy-csv", nargs="+", required=True)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--keep-full-identifiers", action="store_true")
    p.add_argument("--keep-identifier-versions", action="store_true")
    a = p.parse_args()

    tax_assign = MultiLineageDB.load(
        a.taxonomy_csv,
        keep_full_identifiers=a.keep_full_identifiers,
        keep_identifier_versions=a.keep_identifier_versions,
    )
    gather_csvs = tax_utils.collect_gather_csvs(a.gather_csv, from_file=a.from_file)
    results = tax_utils.check_and_load_gather_csvs(
        gather_csvs,
        tax_assign,
        keep_full_identifiers=a.keep_full_identifiers,
        keep_identifier_versions=a.keep_identifier_versions,
    )
    if not results:
        sys.exit("No gather results loaded")
    for r in results:
        r.build_summarized_result()
    with FileOutputCSV(a.output) as fp:
        tax_utils.write_summary(results, fp)
    print(f"{len(results)} queries written to {a.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
