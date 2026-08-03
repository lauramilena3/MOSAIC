#!/bin/bash

# Define the prodigal_multithreaded function
prodigal_multithreaded () {
	    local fasta=$1
	        name="${fasta%.fasta}"  # Remove the .fasta extension
		    original_basename="${fasta%_seq*.fasta}"  # Extract the basename of the original file

		        # Run prodigal on each of the split files
			    prodigal -i ${fasta} -o ${name}_coords_prodigal.txt -a ${name}_ORFs.faa -p meta
		    }

# Export the function so that it can be used by parallel
export -f prodigal_multithreaded

# Get the basename of the input fasta file (without extension)
input_fasta=$1
basename_input="${input_fasta%.fasta}"

# Create a temporary directory for intermediary files
temp_dir="${basename_input}_temp_prodigal"
mkdir -p "$temp_dir"

# Split the input fasta file into smaller fasta files with the desired naming convention
split -l 1000 -d --additional-suffix=.fasta $input_fasta "$temp_dir/${basename_input}_seq"

# Run prodigal_multithreaded function in parallel on the split fasta files
ls $temp_dir/${basename_input}_seq*.fasta | parallel --lb --jobs 64 prodigal_multithreaded {}

# Combine the resulting files into a single file
cat $temp_dir/${basename_input}_seq*_coords_prodigal.txt > ${basename_input}_coords.txt
cat $temp_dir/${basename_input}_seq*_ORFs.faa > ${basename_input}_ORFs.faa

# Remove the intermediary files by removing the entire temporary directory
rm -r "$temp_dir"