export MC=$ILAFORES/analysis/MethodsComparison
cd $MC
source scripts/myFunctions.sh

# Sequence counts by sample
ml seqkit 
:>$MC/Out/stats/all_read_counts_raw.tsv

fastq_files=($(find ./data/*/preproc -maxdepth 2 -type f -name "*paired_1.fastq*"))
total_files=${#fastq_files[@]}

start=0
chunk_size=40

while [ $start -lt $total_files ]; do
    seqkit stats ${fastq_files[@]:$start:$chunk_size} --threads $chunk_size --skip-err >> $MC/Out/stats/all_read_counts_raw.tsv
    start=$((start + chunk_size))
done

grep DNA $MC/Out/stats/all_read_counts_raw.tsv | sponge $MC/Out/stats/all_read_counts_raw.tsv

# Replace first column with dataset and sample IDs :
awk '{
    path = $1;
    rest = substr($0, length(path) + 2);
    split(path, parts, "/");
    group = parts[2];
    sample = parts[4];
    print group " " sample " " rest;
}' OFS='\t' $MC/Out/stats/all_read_counts_raw.tsv > $MC/Out/stats/all_read_counts.tsv

# Summarise mean ± SD by dataset in awk lol just do R gobless deepseek
awk -F' ' '
{
    gsub(/,/, "", $5);  # Remove commas from the fifth column
    group = $1;
    value = $5 + 0;
    
    # Initialize min/max for new groups
    if (!(group in count)) {
        min[group] = value;
        max[group] = value;
    }
    
    # Update min/max
    if (value < min[group]) min[group] = value;
    if (value > max[group]) max[group] = value;
    
    # Track sum and sumsq for mean/stddev
    count[group]++;
    sum[group] += value;
    sumsq[group] += value * value;
}
END {
    for (g in count) {
        mean = sum[g] / count[g];
        stddev = sqrt((sumsq[g] - sum[g]^2 / count[g]) / count[g]);
        printf "%s\t%.2f ± %.2f [%d, %d]\n", g, mean, stddev, min[g], max[g];
    }
}' $MC/Out/stats/all_read_counts.tsv > $MC/Out/stats/sample_count_summary.tsv

