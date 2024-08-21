import sys

def process_fasta(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                # Remove everything after the first space in the header
                header = line.split(" ")[0]
                
                # Split the header at the last underscore
                parts = header.rsplit("_", 1)
                
                # Replace all underscores with hyphens in the first part
                modified_header = parts[0].replace("_", "-")
                
                # Reconstruct the header, keeping the last underscore and part intact
                if len(parts) > 1:
                    modified_header += "_" + parts[1]
                
                # Write the modified header to the output file
                outfile.write(modified_header + "\n")
            else:
                # Write sequence lines unchanged
                outfile.write(line)

if __name__ == "__main__":
    # Check if the correct number of arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python process_fasta.py <input_file> <output_file>")
        sys.exit(1)

    # Get the input and output file names from the command-line arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Process the FASTA file
    process_fasta(input_file, output_file)
