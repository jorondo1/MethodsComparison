# MethodsComparison TODO (untracked: doc/ is in .gitignore)

Started 2026-11-03.

## Results to regenerate

- [ ] **MetaPhlAn per-sample normalisation fixed** (misc_scripts `b86aac6`, pushed 2026-11-03).
  - **The bug:** `parse_MPA(MPA_data = TRUE)` divided by the sum over all samples. Each sample therefore summed to its own classified fraction instead of 100.
  - **To do:** rerun the figures (objects and diversity regenerated 2026-11-04 with `scripts/0_run_all.R`).
  - **Impact measured** on the saved RDS (MetaPhlAn objects only):
    - Alpha diversity is unchanged: every index normalises per sample.
    - Bray-Curtis old vs new correlation: ≥ 0.997 for P19 Saliva/Gut, NAFLD, RA_Gut and PD; 0.979 for Moss; ~0.93 for AD_Skin.
    - PERMANOVA: R² changes in the 3rd digit, and no p-value crosses 0.05.
    - Procrustes: 540 MetaPhlAn pairs; 52 change by > 0.02 and 8 by > 0.05.
      - Moss, MPA_db2025 vs mOTUs4: 0.65 → 0.52.
      - AD_Skin, MPA vs KB90_GTDB: 0.85 → 0.91. Most other AD_Skin MPA pairs rise by 0.02–0.04; MPA vs KB10_GTDB falls from 0.95 to 0.91.
  - [ ] Check any figure or manuscript text about AD_Skin or Moss procrustes agreement involving MetaPhlAn.
- [ ] Classification rates: rerun the mOTUs block (fixed: `awk -F'\t'` and rate = assigned/total), the MetaPhlAn block (100 − Σ species) and the sylph block: `rate_motus`, `rate_mpa`, `rate_sylph` (`scripts/src/classification_rates.sh`, called from `pipeline_2_profiling.sh`).

- [x] Kraken classification rates: RA_Gut has none (kreports deleted). Recovered for the other datasets from git (HANDOFF 2026-11-05); merge the recovered file on the HPC too.

- [x] **sourmash RefSeq counts are doubled** (found 2026-11-03, see `doc/sourmash_refseq_double_counting.txt`). R joins gathers to `SM_RefSeq_20250528_lineages.csv` on the bare accession number, which matches the GCA/GCF twins and versions, so `left_join` duplicates most genomes (+111–116 % on P19_Saliva GQ10–12). Counts, rarefaction depth and count-based alpha are affected. Fixed on the `merged_tables` branch (RefSeq lineages from the sylph-tax table). The relative abundances shift too (max Shannon difference 0.1–0.7 per dataset), so the sourmash RefSeq results need regenerating.

- [ ] **sourmash GTDB 220 missing for Oilfield, Grapevine, Horse** (found 2026-11-10): their objects are `SM_gtdb-rs220-reps`. `parse_merged.R` now renames for every dataset; rerun `0_run_all.R` on the HPC with `rerun <- c('Oilfield', 'Grapevine', 'Horse')`, sync, then rerun the figures.
- [ ] **P19 GQ19 (Saliva) and GQ38 (Gut): the production sourmash sketches come from older read files** (found 2026-11-20, sourmash version check). Production sketches: 2024, from `/nfs3_ib/.../preproc/<s>/<s>_unmatched_2.fastq` (uncompressed). The current `preproc/*.fastq.gz`, which the benchmark sketched, gives other sizes (GQ19: 205,897 against 209,816 hashes; GQ38 sketch 4 % larger in production). The other 6 P19 benchmark samples match exactly. To do: check which reads the other tools used for these samples, and whether more P19 samples are affected (compare sketch sizes, or `sig describe` source and hash count, against a fresh sketch). If so, re-sketch and re-gather them (sourmash only, both databases).
- [ ] PD mOTUs 3: 10 rate rows are NA in `motus_classification_rate.tsv` though the profiles have counts: the rerun of the mOTUs rate block is still pending (item above). Visible as `unknown` in `Out/dataset_summary/dataset_summary.tsv`.

## Methods notes (for the manuscript)

- [ ] **Adapter clipping was skipped in 9 datasets (found 2026-11-18). Decide what to do (the user's call).**
  - In the HPC kneaddata runs since 2024, `ILLUMINACLIP` pointed at `/cvmfs/soft.mugqic/.../TruSeq3-PE-2.fa`. The kneaddata container mounts only the roots of the host-index paths (`/net`), so Trimmomatic couldn't read the file. It clipped no adapters, raised no error, and still applied `SLIDINGWINDOW:4:20 MINLEN:50`.
  - Test: a kneaddata log with adapters loaded has 4 `Using Long Clipping Sequence` lines. Counts: P19_Gut 0/28, NAFLD 0/35, AD_Skin 0/148, PD 0/724, Bee 0/40, RA_Gut 0/168, Maple 0/88, Oilfield 0/87, Grapevine 0/78. Adapters were clipped in P19_Saliva 25/26 and Moss 140/140 (2023 runs) and Horse 744/744 (Fir, adapter file copied under /scratch).
  - Size of the effect, on Horse SRR28771685 (same raw reads, with vs without adapters): 0.3 % more pairs survive trimming without adapter clipping (547,965 vs 546,276). Read-through adapter tails stay on the reads.
  - The methods currently say Trimmomatic with adapter clipping for every dataset, so either the text changes or the 9 datasets are re-QC'd and reprofiled.
  - Fix for any new run: point `--trimmomatic_adapters` at `$FAST/adapters/TruSeq3-PE-2.fa`, a copy of the /cvmfs file (same md5) under `/net`, which the container mounts. Check one log for the 4 `Using Long Clipping Sequence` lines.

- [ ] **Response to reviewers: sample numbers (a reviewer found them incoherent across the ms).** Two sample sets, stated once in the methods and in the dataset table (`Out/dataset_summary/dataset_table.html`, `scripts/6.2_dataset_summary.R`):
  - every post-QC sample for the comparisons that need no metadata (alpha diversity, PCoA, distances, Procrustes; 2,216 samples);
  - only samples with a value for the grouping variable for the Wilcoxon tests and PERMANOVA (2,185): Moss 131 of 140, Maple 66 of 88.
  - Removed from everything beforehand (`scripts/src/dropped_samples.tsv`):
    - Gut COV: 7 controls (mock community MSA-3001, GQ27–29, GQ54–56);
    - PD: 5 runs that are neither PD nor control (library SP*/SC*);
    - Skin AD: 78 runs absent from the study's metadata, mostly very shallow (only the 70 with metadata are used).
  - Figure labels must use the matching n.
  - [ ] Takes effect in the results after the HPC rerun of Moss (HANDOFF).
- [ ] **The sample sets differ between the diversity comparisons and the statistical tests.**
  - The diversity comparisons use every sample of a dataset: alpha diversity, PCoA, pairwise distances and Procrustes.
  - The Wilcoxon tests and PERMANOVA use only the samples that have a value for the grouping variable (`grouping_variable` in `0_config.R`).
  - Maple: 88 samples in the diversity comparisons, but only the 66 with a Provenance (33 M, 33 S) in the tests.
  - The `n` in the figure labels (`dataset_names_n_lookup_alphadiv.RDS`) comes from the Wilcoxon tests, so for Maple it's the tested count.
  - State this in the methods.
- [ ] **PD mOTUs 3: 70 profiles salvaged (tell reviewers).** 70 of the 725 PD mOTUs 3 profiles (rows 1–25 plus 45 scattered rows) were empty files from a batch that died on 2025-03-29, so the merged table had 655 samples. Until 2026-11-07 the classification-rate script also gave each of them the previous sample's rate. Rerun on 2026-11-07 (job 407787, same mOTUs 3.1.0 and database). Mention the rerun in the methods or the response to reviewers, then rerun `merge_tables.sh PD` and the PD mOTUs 3 figures.
- [ ] **mOTUs left out of the classification rates (decided 2026-11-09): say why in the methods or the ms.** mOTUs only sees reads of 10 single-copy prokaryotic marker genes. Its assigned/(assigned + unassigned) is therefore a share of the prokaryotic marker-gene signal, not of the reads, and it ignores host, eukaryotic and non-marker DNA. It is not comparable with the other tools (Bee median 0.95, against 0.11 for MetaPhlAn).
- [ ] **Kraken identifiability = species-rank share (2026-11-09): say so in the methods.** The rate is the fraction of reads Kraken itself assigns at species rank or below, read from Bracken's summary. It is not the Bracken total, which redistributes nearly all classified reads to species, and not the any-rank classified rate the figure used before. Example, PD at confidence 0.10 (GTDB): species-rank 0.55, any rank 0.99. At 0.90 it is ~0 everywhere. The rates are in `Out/classification_rates/kraken_species_classification_rate.tsv` (`rate_kraken_species`). Caveat: it is per read, so reads from conserved regions of known species stay at higher ranks.
  - [x] Mac: Kraken switched to this file in `scripts/6.1_Identifiability.R` (2026-11-10; replaces `8_Classification_rates.R`). It has the same 4 columns and path layout as `kraken_classification_rate.tsv`.
- [ ] **Discussion: at Kraken confidence 0.45, RefSeq puts a larger share of reads at species rank than GTDB (human datasets), the reverse of 0.10 and of every other tool** (checked 2026-11-10; not a parsing bug: the Bracken-summary totals equal the kreport totals for all 13,819 samples).
  - **The numbers** (medians, human datasets, 0.45): GTDB classifies 0.87–0.93 of the reads at any rank but only 0.07–0.16 at species; RefSeq classifies 0.49–0.72 at any rank and 0.13–0.25 at species. Within the classified reads, the species share is 0.08–0.18 for GTDB against 0.21–0.53 for RefSeq. In non-human datasets (Horse, soils, Maple) GTDB is higher at every confidence, because RefSeq barely classifies anything there.
  - **Likely explanation:** Kraken assigns a read to the deepest node supported by ≥ confidence × its k-mers. GTDB r220 has ~113k species with one representative genome each and splits gut genera finely (e.g. *Bacteroides*, *Prevotella*). A read's k-mers are spread over sister species, or moved up to the genus by the LCA step, so species-level support rarely reaches 0.45 and the read stops at genus. RefSeq has fewer species but many strains of the common human-associated ones, so a read either clearly matches a well-covered species or stays unclassified. A denser database classifies more reads but at higher ranks: Nasko et al. 2018, Genome Biol 19:165 (RefSeq growth and k-mer LCA species identification).
  - **For the ms:** the species-rank share mixes "is the read classified?" with "how deep?". It is the proxy for "share of reads identified at species level", but higher confidence trades one against the other differently in each database.
  - [ ] Optional HPC check: recompute the species-rank reads from one kreport per database (Grapevine KB45 vs KB45_GTDB; awk one-liner in the 2026-11-10 session) and compare with the rate file.
- Out of scope for this project: rarefaction curves (`9.1`/`9.2`) and the DAA scripts (`5.0`–`5.4`).
- [ ] **New datasets (Nov 2026): describe them in the methods.**
  - **Shared preprocessing:** raw paired reads were retrieved from the NCBI Sequence Read Archive. Each run's read count was checked against the archive's.
    - QC used KneadData 0.12.0, with Trimmomatic (`SLIDINGWINDOW:4:20 MINLEN:50`), then removal of host reads by mapping to the host genome with Bowtie2 (`--very-sensitive-local`).
    - These settings are the same as for Maple.
  - **Oilfield (PRJNA1196024):** oil-contaminated coastal wetland soils, surface layer (0–20 cm).
    - Paper: Commun. Earth Environ. 2025, doi:10.1038/s43247-025-02062-1.
    - 87 samples, one per site: Shengli Oilfield 44, Dagang Oilfield 43. Illumina NovaSeq 6000, 2×150.
    - **Grouping variable:** oilfield, taken from the sample names (`-SL` / `-DG`).
    - **Host removal:** human GRCh38 (contamination filter only; no host).
  - **Grapevine (PRJNA1153634):** root rot of wine grapevines, Ningxia, China.
    - Paper: Genome Biol. 2026, doi:10.1186/s13059-025-03905-y.
    - 78 samples from 39 vines: healthy vs diseased × 2 sites × 2 cultivars (Cabernet Sauvignon, Syrah) × rhizosphere soil and root. DNBSEQ-T7, 2×150.
    - Metadata decoded from the sample names (e.g. `UHDZR3`).
    - **Grouping variable:** health status or compartment (to decide).
    - **Pairing:** each vine contributes a soil and a root sample, so the compartment comparison is paired. Our tests are unpaired, so state it as a limitation, or test one compartment.
    - **Host removal:** *Vitis vinifera* PN40024 T2T reference (GCF_030704535.1).
  - **Snake (PRJNA1353898):** cloacal swabs of Eastern Indigo Snakes (*Drymarchon couperi*), Florida, 2025.
    - Paper: PLOS ONE 2026, doi:10.1371/journal.pone.0350824.
    - 80 snakes: 60 captive (one facility) and 20 wild. Illumina NovaSeq X.
    - **Grouping variable:** wild vs captive, from the BioSample location.
    - The paper splits the captive snakes by *Cryptosporidium serpentis* status (30/30), but that status isn't public.
    - **Host removal:** *Drymarchon corais* (GCA_043091225.1), the reference the authors used.
    - Status: only 3 samples are downloaded so far, as a test.

## Adding a dataset

`scripts/0_run_all.R` builds the phyloseq objects and the diversity metrics one dataset at a time. Results are cached per dataset in `Out/_Rdata/cache/<DS>/`, and the script then assembles the usual `Out/_Rdata/*.RDS` files. Only datasets without a complete cache are computed; list datasets in `rerun` to recompute them anyway.

1. HPC:
   - run the profilers;
   - add a line to `scripts/src/datasets.tsv`;
   - run `bash scripts/src/merge_tables.sh <DS>`, plus `scripts/utils/sylph_merge.sh <DS> <db> <sample tsv>` for each sylph database.
2. `scripts/utils/sample_data.R`: add the sample data to `load_sample_data()`'s list.
3. `scripts/0_config.R`: add `grouping_variable` and `dataset_names` entries.
4. Add `<DS>` to `datasets` in `scripts/0_run_all.R` and run it on the HPC (`Rscript scripts/0_run_all.R`).
   - Before parsing, it checks steps 1–3 and reports anything missing, with the fix.
5. Sync (`scripts/utils/sync_from_hpc.sh`), then make the figures.

## sylph

- [ ] Rerun sylph 1.0.0 (env `$ILAFORES/programs/conda_envs/sylph`) for all datasets and DBs. First delete the old `per_sample/` dirs and the sketches made with 0.9.0.
- [ ] Merge with `sylph_merge.sh [-f]`.
- [ ] One-off: build `merged_species.tsv` from the existing `merged_all_ranks.tsv`, then delete `merged_all_ranks.tsv`.
- [ ] Check with the sylph authors that dividing the counts by the classification rate is correct (`doc/sylph_methods_note.txt`).
- [ ] Figures that include sylph: they need a rethink (only the parsing and config are done).

## Code issues noticed

- [x] `1_diversity_metrics.R` sourcing `myFunctions.R` / `0_Config.R`: fixed. It runs on the HPC.
- [x] `sync_from_hpc.sh`: now tracked, as `scripts/utils/sync_from_hpc.sh`.

## Merged tables (done, branch `merged_tables`)

All tools are read from the merged tables (`scripts/src/merge_tables.sh`, `scripts/utils/parse_merged.R`). The objects were checked against the per-sample parsing before the comparison scripts were archived. They were identical, apart from MetaPhlAn (float precision) and sourmash RefSeq (the double-counting fix).

## Before the next Zenodo snapshot

- [ ] Purge data/ and Out/ from the git history: `doc/git_history_purge_steps.txt`.

## To launch

- [x] **Horse on Fir** (PRJNA1102860): profiled and merged 2026-11-07, outputs in `data/Horse` on Fir (2026-11-09). Remaining: the R steps of "Adding a dataset". Original note: Merge the `merged_tables` and ILL `kneaddata_gzip` PRs first. Set `SRC` in `setup_fir.sh`, and `HORSE_FASTA` in `horse.sh`, or `HORSE_INDEX` in `setup_fir.sh` if the index is built on iv12. Steps are in `doc/HANDOFF_FIR.md`.
- [x] **Resource benchmark** for the reviewer's runtime and memory request: 55 samples x 17 tools done on cp0359 (2026-11-16). Scripts: `scripts/benchmark/`, figure `scripts/10_benchmark_plot.R`, methods `scripts/benchmark/benchmark_methods.txt`.
  - [ ] Add Horse (rows 56-60) and rerun sourmash with 4.9.4 + branchwater, 24 threads for both databases (RefSeq SBT -> zip). Steps in `doc/HANDOFF.md`.
  - [ ] Sub-agent check: does 4.9.4 change the taxonomic composition (HANDOFF step 6)? If it does, decide whether the main sourmash profiles need a rerun.
  - [ ] Methods text: node hardware, Horse, 60 samples, sourmash 4.9.4.
