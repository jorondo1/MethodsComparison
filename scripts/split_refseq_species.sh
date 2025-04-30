#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH -D /net/nfs-ip34/home/def-ilafores/analysis/MethodsComparison
#SBATCH -o /net/nfs-ip34/home/def-ilafores/analysis/MethodsComparison/logs/split_refseq-%A_%a.slurm.out
#SBATCH --time=2:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=1G
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
echo "Compressing and copying back to ${OUT_DIR}..."

mkdir -p $temp_dir/out
#for species_file in "$temp_dir"/*.fna; do
#    species_name=$(basename "$species_file" .fna)
#    if [[ -s "${OUT_DIR}/signatures/${base_name}_${species_name}.fna.gz" ]]; then
#    continue
#    else
#    gzip -c "$species_file" > "${temp_dir}/out/${base_name}_${species_name}.fna.gz"
#    fi
#done

# Sourmash signatures
find $temp_dir -name '*.fna' > $temp_dir/file_list.txt

echo "Computing sourmash signatures..."
ml apptainer

mkdir -p "$OUT_DIR"/signatures/tmp

singularity exec --writable-tmpfs -e -B $ANCHOR$ILAFORES:$ANCHOR$ILAFORES,$ANCHOR/fast2/def-ilafores:$ANCHOR/fast2/def-ilafores $ANCHOR$ILL_PIPELINES/containers/sourmash.4.8.11.sif bash -c "
for FNA_file in $(find $temp_dir -name '*.fna'); do
    # sketch signature
    sourmash sketch dna -p k=31,scaled=1000,abund --name-from-first $FNA_file --outdir $temp_dir/out/
    
    # check signature
    sourmash signature describe $temp_dir/out/${FNA_file}.sig
    
    # Copy signature
    cp $temp_dir/out/${FNA_file}.sig \"$OUT_DIR\"/signatures/tmp
    mv \"$OUT_DIR\"/signatures/tmp/${FNA_file}.sig \"$OUT_DIR\"/signatures
    rm $temp_dir/out/${FNA_file}.sig
done
"

#mkdir -p "$OUT_DIR"/genomes

# cp $temp_dir/out/*fna.gz "$OUT_DIR"/genomes

# Clean up
rm -rf "$temp_dir"
rm -r "$OUT_DIR"/signatures/tmp
echo "Done!"
