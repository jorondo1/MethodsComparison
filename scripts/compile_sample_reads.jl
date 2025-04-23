#!/usr/bin/env julia

# Compile the number of reads in metagenome samples

# ========== PACKAGES ==========
using ArgParse, Distributed, Glob
using FASTX, DataFrames, CSV 

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
			help = "Only count paired samples (will look for the string 'paired' in filenames)"
			default = true
		"input_directories"
			help = "Input directories to process (at least one required)"
			required = true  
			arg_type = String
			nargs = '+'  # "+" means one or more arguments
	end
	return parse_args(s)
end

# ========== PARALLEL FUNCTION SETUP ==========
function init_workers(ncores)
    addprocs(ncores)
    @everywhere begin
        using FASTX
        # ========== COUNT READS ==========
        function count_reads_fastx(filename)
            try
                reader = FASTQ.Reader(open(filename))
                count = 0
                for _ in reader
                    count += 1
                end
                close(reader)
                return count
            catch e
                @warn "Failed to process $filename: $e"
                return missing
            end
        end
    end
end

# ========== GENERATE FILE LIST ==========
function generateFastaList(directories::Vector{String})
    fasta_files = String[]
    for dir in directories
        files = Glob.glob("*/*.fa*", dir)
        for file in files # only keep 1st mate in pair, without contam files
            if occursin("_paired_", file) &&
               !(occursin("contam", file) || occursin("_2.", file) || occursin("unmatched", file))
                push!(fasta_files, file)
            end
        end
    end
    return(fasta_files)
end

# ========== GENERATE REPORT ==========
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

# ========== MAIN ==========
function main()
    args = parse_commandline()
    init_workers(args["ncores"])
    output_file_path = joinpath(args["output_dir"], "read_counts.csv")
    generate_read_counts(args["input_directories"], output_file_path)
end

main()