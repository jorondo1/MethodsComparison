#!/usr/bin/env julia

# =============================================================================
# compile_sample_reads.jl
#
# Counts reads in FASTQ(.gz) files across multiple input directories,
# then produces:
#   1. <output_dir>/all_read_counts.tsv     — per-sample counts
#   2. <output_dir>/sample_count_summary.tsv — mean ± SD [min, max] per dataset
#   3. <output_dir>/failed_samples.tsv      — files that could not be parsed
#
# Parallels the logic of count_reads.sh but runs entirely in Julia,
# using worker processes for fast parallel I/O.
# =============================================================================

# ========== PACKAGES ==========
using ArgParse, Distributed, Glob
using DataFrames, CSV, Statistics

# ========== ARG PARSING ==========
function parse_commandline()
    s = ArgParseSettings(
        description = "Count reads in FASTQ samples and summarise by dataset."
    )

    @add_arg_table s begin
        "--ncores", "-c"
            help    = "Number of parallel worker processes"
            arg_type = Int
            default  = 4
        "--output_dir", "-o"
            help     = "Directory where output TSVs will be written"
            arg_type = String
            required = true
        "--paired_only"
            help    = "Only count files whose name contains 'paired' (flag, no value needed)"
            action  = :store_true   # present → true, absent → false
        "--glob_pattern"
            help    = "Glob pattern used to find FASTQ files inside each input directory"
            arg_type = String
            default  = "*/*.fastq*"   # matches one sub-level deep; adjust if needed
        "input_directories"
            help     = "One or more root directories to search for FASTQ files"
            required = true
            arg_type = String
            nargs    = '+'
    end

    return parse_args(s)
end

# ========== PARALLEL SETUP ==========
# Parse args first so we know how many workers to spin up.
args = parse_commandline()

# addprocs launches extra Julia processes; each one will load the packages below.
addprocs(args["ncores"])

# @everywhere runs the block on the main process AND every worker.
@everywhere using FASTX, CodecZlib

# ========== READ-COUNTING FUNCTION ==========
# This function runs on whichever worker picks up the file.
# It returns the integer count on success, or `missing` on failure.
@everywhere function count_reads_fastq(filename::String)
    try
        reader = if endswith(filename, ".gz")
            # Wrap the file in a gzip decompressor before handing to FASTQ.Reader
            FASTQ.Reader(GzipDecompressorStream(open(filename)))
        else
            FASTQ.Reader(open(filename))
        end

        # Iterate through every record and count — never loads all records into RAM
        n = 0
        for _ in reader
            n += 1
        end
        close(reader)
        return n

    catch e
        # Return the error message as a string so the caller can log it
        return "ERROR: $(e)"
    end
end

# ========== FILE DISCOVERY ==========
# Finds FASTQ files under each input directory using a glob pattern,
# then optionally filters to paired-end first mates only.
function find_fastq_files(directories::Vector{String},
                          glob_pattern::String,
                          paired_only::Bool)

    found = String[]
    for dir in directories
        # Glob.glob searches relative to `dir`
        for f in Glob.glob(glob_pattern, dir)
            if paired_only
                # Keep only first-mate files; skip second mate, unmatched, contam
                keep = occursin("paired", f) &&
                       !occursin("_2.", f)   &&
                       !occursin("unmatched", f) &&
                       !occursin("contam", f)
                keep || continue
            end
            push!(found, f)
        end
    end

    isempty(found) && @warn "No FASTQ files found — check your input directories and --glob_pattern."
    return found
end

# ========== PATH → (dataset, sample) PARSING ==========
# Mirrors the awk in count_reads.sh:
#   ./data/DatasetName/preproc/SampleID/SampleID_paired_1.fastq.gz
#          parts[2]                parts[4]
# Adjust the indices below if your directory layout differs.
function parse_path(filepath::String)
    parts = splitpath(filepath)
    # parts[1] is likely "." or an absolute root; data is typically at parts[2]
    dataset = length(parts) >= 2 ? parts[2] : "unknown"
    sample  = length(parts) >= 4 ? parts[4] : basename(filepath)
    return dataset, sample
end

# ========== MAIN PIPELINE ==========
function main()
    # ------------------------------------------------------------------
    # 1. Discover files
    # ------------------------------------------------------------------
    fastq_files = find_fastq_files(
        args["input_directories"],
        args["glob_pattern"],
        args["paired_only"]
    )

    n_files = length(fastq_files)
    println("Found $n_files FASTQ file(s). Counting reads on $(args["ncores"]) workers…")

    # ------------------------------------------------------------------
    # 2. Count reads in parallel
    #    @distributed (vcat) merges results from all workers into one vector.
    # ------------------------------------------------------------------
    results = @distributed (vcat) for file in fastq_files
        [(filepath = file, result = count_reads_fastq(file))]
    end

    # ------------------------------------------------------------------
    # 3. Separate good results from failures
    # ------------------------------------------------------------------
    good_rows  = NamedTuple{(:dataset, :sample, :filepath, :num_seqs), Tuple{String,String,String,Int}}[]
    error_rows = NamedTuple{(:filepath, :error), Tuple{String,String}}[]

    for (; filepath, result) in results
        if result isa Int
            dataset, sample = parse_path(filepath)
            push!(good_rows, (dataset = dataset,
                              sample  = sample,
                              filepath = filepath,
                              num_seqs = result))
        else
            # result is an error string
            push!(error_rows, (filepath = filepath, error = string(result)))
        end
    end

    # ------------------------------------------------------------------
    # 4. Write all_read_counts.tsv
    # ------------------------------------------------------------------
    mkpath(args["output_dir"])
    counts_path = joinpath(args["output_dir"], "all_read_counts.tsv")

    counts_df = DataFrame(good_rows)
    sort!(counts_df, [:dataset, :sample])
    CSV.write(counts_path, counts_df; delim = '\t')
    println("  → Wrote $(nrow(counts_df)) rows to $counts_path")

    # ------------------------------------------------------------------
    # 5. Write sample_count_summary.tsv  (mean ± SD [min, max] per dataset)
    # ------------------------------------------------------------------
    summary_path = joinpath(args["output_dir"], "sample_count_summary.tsv")

    summary_rows = []
    for (dataset, group_df) in pairs(groupby(counts_df, :dataset))
        vals = Float64.(group_df.num_seqs)
        m   = mean(vals)
        sd  = std(vals; corrected = false)   # population SD, matching the awk script
        mn  = minimum(vals)
        mx  = maximum(vals)
        push!(summary_rows, (
            dataset = dataset.dataset,
            summary = @sprintf("%.2f ± %.2f [%d, %d]", m, sd, Int(mn), Int(mx))
        ))
    end

    # @sprintf needs Printf
    using Printf   # loaded lazily here so it's only needed when we reach this point
    summary_df = DataFrame(summary_rows)
    sort!(summary_df, :dataset)
    CSV.write(summary_path, summary_df; delim = '\t')
    println("  → Wrote $(nrow(summary_df)) dataset summaries to $summary_path")

    # ------------------------------------------------------------------
    # 6. Write failed_samples.tsv (if any)
    # ------------------------------------------------------------------
    if !isempty(error_rows)
        errors_path = joinpath(args["output_dir"], "failed_samples.tsv")
        errors_df   = DataFrame(error_rows)
        CSV.write(errors_path, errors_df; delim = '\t')
        @warn "$(nrow(errors_df)) file(s) failed — see $errors_path"
    else
        println("  → No failures detected.")
    end

    println("Done.")
end

main()
