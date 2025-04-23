#!/usr/bin/env julia

# Compile the number of reads in metagenome samples

#= #............
Parse arguments.
=# #............
using ArgParse

function parse_commandline()
	s = ArgParseSettings()
    
	@add_arg_table s begin
		"--n_cores", "-c"
			help = "Number of cores for parallel processing"
			arg_type = Int
			default = 2
		"--output_dir", "-o"
			help = "Path where reports will be saved"
			arg_type = String
#			required = true
		"--paired_only"
			help = "Only count paired samples (will look for the string 'paired' in filenames)"
			default = true
		"input_directories"
			help = "Input directories to process (at least one required)"
#			required = true  
			arg_type = String
			nargs = '+'  # "+" means one or more arguments
	end
	
	return parse_args(s)
	
end

#= #...................................
Adding workers for parallel processing.
=# #...................................

#= *************************
Function Generate file list*
=# #************************
using Glob
function generateFastaList(directories::Vector{String})
    fasta_files = String[]
    for dir in directories
        files = Glob.glob("*/*.fa*", dir)
        for file in files
            # Keep ONLY if:
            # 1. Contains "_paired_" AND
            # 2. Does NOT contain "contam", "_2.", or "unmatched"
            if occursin("_paired_", file) &&
               !(occursin("contam", file) || occursin("_2.", file) || occursin("unmatched", file))
                push!(fasta_files, file)
            end
        end
    end
    println(fasta_files)
end

#= ******************
Function Count reads*
=# #*****************

#= **********************
Function Generate report*
=# #*********************

#= 
Main 
=# 

function main()
	args = parse_commandline()
	generateFastaList(args["input_directories"])
end

main()