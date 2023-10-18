#!/bin/bash

# Define the prodigal_multithreaded function
prodigal_multithreaded () {
        local fasta=$1
        name="${fasta%.fasta}"
        echo $name
        prodigal -i ${fasta} -o ${name}_coords_prodigal.txt -a ${name}_ORFs.faa -p meta
}

# Export the function so that it can be used by parallel
export -f prodigal_multithreaded

# Split the input fasta file into smaller fasta files

split -l 2000 -d --additional-suffix=.fasta $1 seq

# Run prodigal_multithreaded function in parallel on the split fasta files
ls seq*fasta | parallel --lb --jobs 64 prodigal_multithreaded {}