#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH -D /net/nfs-ip34/home/def-ilafores/analysis/MethodsComparison
#SBATCH -o /net/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/logs/split_refseq-%A_%a.slurm.out
#SBATCH --time=2:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=24G
#SBATCH -A def-ilafores
#SBATCH -J split_refseq

export ANCHOR="/net/nfs-ip34"
export FNA_PATH=${ANCHOR}$(awk "NR==$SLURM_ARRAY_TASK_ID" "${ANCHOR}${1}")
export OUT_DIR="${ANCHOR}${2}"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Get the base filename without extension
base_name=$(basename "$FNA_PATH" .fna.gz)

# Create a temporary directory for this file
temp_dir=$(mktemp -d -p /tmp "split_${SLURM_JOB_ID}_XXXXXX") || exit 1

# Stagger jobs
echo "Copying $FNA_PATH to ${temp_dir}... "
sleep $((RANDOM % 30))
cp "$FNA_PATH" $temp_dir

echo "Splitting file..."
zcat $temp_dir/$(basename "$FNA_PATH") | awk -v temp_dir="$temp_dir" '
/^>/ {
    # Extract the species identifier (first 9 characters after >)
    species_id = substr($1, 1, 10)
    
    # Remove the > character
    species_id = substr(species_id, 2)
    
    # Create filename
    output_file = temp_dir "/" species_id ".fna"
    
    # Print to the appropriate file
    print > output_file
    next
}
{
    # Print sequence lines to the same file
    print >> output_file
}'

# Compress and move the output files

echo "Compress and copy back to ${OUT_DIR}..."
for species_file in "$temp_dir"/*.fna; do
    species_name=$(basename "$species_file" .fna)
    gzip -c "$species_file" > "${OUT_DIR}/${base_name}_${species_name}.fna.gz"
done

# Clean up
rm -rf "$temp_dir"
echo "Done!"
