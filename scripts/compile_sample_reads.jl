#!/usr/bin/env julia

# ========== TOP LEVEL PACKAGES (Main process only) ==========
using ArgParse, Glob, DataFrames, CSV

# ========== ARG PARSING ==========
function parse_commandline()
    s = ArgParseSettings()
    
    @add_arg_table s begin
        "--ncores", "-c"
            help = "Number of cores for parallel processing"
            arg_type = Int
            default = 2
        "--output_dir", "-o"
            help = "Path where reports will be saved"
            arg_type = String
            required = true
            default = "."
        "--paired_only"
            help = "Only count paired samples (will look for 'paired' in filenames)"
            default = true
        "input_directories"
            help = "Input directories to process (at least one required)"
            required = true  
            arg_type = String
            nargs = '+'
    end
    return parse_args(s)
end

# ========== PARALLEL WORKER SETUP ==========
function setup_workers(ncores)
    using Distributed
    addprocs(ncores)
    
    # Load required packages on all workers
    @everywhere using FASTX
    
    # Define counting function on all workers
    @everywhere function count_reads_fastx(filename)
        try
            FASTQ.Reader(open(filename)) do reader
                sum(1 for _ in reader)
            end
        catch e
            @warn "Failed: $filename ($e)"
            missing
        end
    end
end

# ========== FILE PROCESSING FUNCTIONS ==========
function generateFastaList(directories::Vector{String})
    fasta_files = String[]
    for dir in directories
        files = Glob.glob("*/*.fa*", dir)
        for file in files
            if occursin("_paired_", file) &&
               !(occursin("contam", file) || occursin("_2.", file) || occursin("unmatched", file))
                push!(fasta_files, file)
            end
        end
    end
    return fasta_files
end

function generate_read_counts(input_dirs::Vector{String}, output_path::String)
    fastq_files = generateFastaList(input_dirs)  
    
    # Parallel counting
    results = @distributed (vcat) for file in fastq_files
        (filepath = file, count = count_reads_fastx(file))
    end

    # Save as CSV
    df = DataFrame(results)
    CSV.write(output_path, df)
    return df
end

# ========== MAIN EXECUTION ==========
function main()
    args = parse_commandline()
    setup_workers(args["ncores"])
    output_file_path = joinpath(args["output_dir"], "read_counts.csv")
    generate_read_counts(args["input_directories"], output_file_path)
end

main()