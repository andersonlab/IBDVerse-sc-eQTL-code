#!/bin/bash
module load HGI/softpack/groups/team152/genetics_tools/1

# Define the window size
window=2e6

# Directory where the files are located
file_dir="." # Change this if files are in a different directory
output_dir="$file_dir/collapsed" # Change this if you want to output to a different directory
mkdir -p "$output_dir"

# Use find to list all files in the current directory and extract unique condition patterns
arr=($(ls -p $file_dir | grep -v / | cut -f 4 -d '.' | sort -u))
echo "Conditions found: ${arr[@]}"

# Extract unique trait patterns (like IBD, UC)
arr2=($(ls -p $file_dir | grep -v / | cut -f 1 -d '.' | sort -u))
echo "Traits found: ${arr2[@]}"

# Loop through conditions and traits and concatenate files
for cond in "${arr[@]}"
do
   for trait in "${arr2[@]}"
   do
      # Initialize a temporary file for concatenation
      tmp_file="${output_dir}/${trait}.featureCounts.${window}.${cond}.txt"
      > "$tmp_file"

      # Check each chunk file before concatenating
      for chunk in {1..20}
      do
         chunk_file="${file_dir}/${trait}.featureCounts.${window}.${cond}.chunk.${chunk}.txt"
         if [[ -f "$chunk_file" ]]; then
            awk 'NR>1' "$chunk_file" >> "$tmp_file" # Skip header when concatenating
         else
            echo "Skipping non-existent file: $chunk_file"
         fi
      done
   done
done

# Tar, compress, and remove individual chunk files
for trait in "${arr2[@]}"
do
  tar_cmd="tar -zcvf ${output_dir}/${trait}.tgz ${output_dir}/${trait}.featureCounts.${window}.*.chunk.*.txt --remove-files"
  echo $tar_cmd
  eval $tar_cmd
  gzip_cmd="cat ../coloc_header.txt ${output_dir}/${trait}.featureCounts.${window}* | bgzip -c > ${output_dir}/${trait}.gz"
  echo $gzip_cmd
  eval $gzip_cmd
  rm_cmd="rm ${output_dir}/${trait}.featureCounts.${window}*"
  eval $rm_cmd
done
