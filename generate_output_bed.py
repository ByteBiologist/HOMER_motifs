import gzip
import subprocess
import os
import concurrent.futures
import sys

genome_build = sys.argv[1]

input_directory = f"/mnt/ebs/jackal/FILER2/FILER2-production/Homer/motif/split_by_motifName/{genome_build}/"
output_directory = f"/mnt/ebs/jackal/FILER2/FILER2-production/Homer/motif/processed_bed/{genome_build}/"
#output_directory = f"/mnt/ebs/jackal/FILER2/FILER2-production/Homer/motif/"

def get_reference(genome_build, sequenceID):
    if genome_build == "hg19":
        fa = "hg19.fa"
    elif genome_build == "hg38":
        fa = "hg38.fa"
    else:
        raise ValueError("Invalid genome build specified")

    output = subprocess.check_output(["samtools", "faidx", fa, sequenceID]).decode("utf-8")
    lines = output.split('\n')
    ref_sequence = ''.join(lines[1:])
    ref_sequence = ref_sequence.lower() # Convert to lowercase
    return ref_sequence

def get_cellandmotif(name, length, column_mapping):
    key_with_length = (name, length)  # Create the key using name and length

    if key_with_length in column_mapping:
        key = key_with_length
    else:
        key = (name, None)  # Create a key without length to prioritize name match

    if key in column_mapping:
        motif_id = column_mapping[key][0]
        cell_type = column_mapping[key][1]
        DNA_binding = column_mapping[key][2]
        consensus_length = column_mapping[key][3]
    else:
        motif_id = "."
        cell_type = "."
        DNA_binding = "."
        consensus_length = "."

    return motif_id, cell_type, DNA_binding, consensus_length

def get_reverse_complement(sequence):
    complement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    reverse_complement = ''
    sequence = sequence.lower()  # Convert to lowercase
    for base in reversed(sequence):
        if base in complement:
            reverse_complement += complement[base]
        else:
            reverse_complement += base
    return reverse_complement

# Load the data dictionary from parsed_subheadings.txt
column_mapping = {}

with open("parsed_subheadings.txt", "r") as parsed_file:
    next(parsed_file)  # Skip the header line
    for line in parsed_file:
        columns = line.strip().split("\t")
        # Extract the keys (columns[1] and columns[12]) and value (columns 0 and 4) from the current line
        key = (columns[1], columns[10])
        value = [columns[0], columns[4], columns[3], columns[9]]
        # Add the key-value pair to the dictionary
        column_mapping[key] = value

# Create the output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

def process_file(input_filepath):
    input_filename = os.path.basename(input_filepath)  # Extract the input filename
    output_filename = input_filename.replace(".bed.gz", ".processed.bed")
    output_filepath = os.path.join(output_directory, output_filename)

    # Check if the output file (with .gz extension) exists in the output directory
    if os.path.exists(output_filepath + ".gz"):
        print(f"Output file {output_filepath}.gz already exists. Skipping.")
        return

    #  Check if the current file matches the three specified filenames
    if input_filename in ["FOXA1(Forkhead).bed.gz", "Nr5a2(NR).bed.gz", "OCT:OCT(POU,Homeobox).bed.gz"]:
        print(f"Running motif_score.py for {input_filename}")
        subprocess.run(["python3", "motif_score.py", input_filepath, output_directory, genome_build])
        return

    # Open the output file for writing
    with open(output_filepath, "w", encoding="utf-8") as output_file:
#        output_file.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tmotif_length\tmotifID\tDNA_binding_domain\tcell_type\tbinding_sequence\tconsensus_length\treverse_strand\n")
        output_file.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tmotifID\tDNA_binding_domain\tcell_type\tbinding_sequence\tconsensus_sequence\n")
        # Open the bed file and process lines on-the-fly
        with gzip.open(input_filepath, "rt", encoding="utf-8") as bed_file:
            for line in bed_file:
                if line.startswith("#"):  # Skip comment lines
                    continue

                columns = line.strip().split("\t")
                motif_name = columns[3]  # Extract motif name

                start = int(columns[1])
                end = int(columns[2])
                motif_length = end - start + 1  # Calculate motif_length, 1-based coordinate system
                strand = columns[5]
                sequenceID = f"{columns[0]}:{columns[1]}-{columns[2]}"

                key = (columns[3], str(motif_length))

                if key in column_mapping:
                    consensus_sequence = column_mapping[key][3]
                    if motif_length != len(consensus_sequence):
                        print(motif_length, "!=", consensus_sequence)
                        # Skip this line as motif_length doesn't match consensus_sequence length
                        continue

                    motif_id, cell_type, DNA_binding, consensus_sequence = get_cellandmotif(columns[3], str(motif_length), column_mapping)

                    if strand == '+':
                        binding_sequence = get_reference(genome_build, sequenceID)
                    elif strand == '-':
                        ref_sequence = get_reference(genome_build, sequenceID)
                        binding_sequence = get_reverse_complement(ref_sequence)

                    output_line = "\t".join([columns[0], str(int(columns[1]) - 1), columns[2], columns[3], columns[4], columns[5]] + [ motif_id, DNA_binding, cell_type, binding_sequence, consensus_sequence])
                    output_file.write(output_line + "\n")

        # Compress the output file using bgzip
        subprocess.run(["bgzip", output_filepath])

        print(f"Output written to {output_filepath}")

        print("Output compressed to", output_filepath + ".gz")

# Process each input file using multiple threads
with concurrent.futures.ThreadPoolExecutor() as executor:
    input_files = [os.path.join(input_directory, filename) for filename in os.listdir(input_directory) if filename.endswith(".bed.gz")]
    executor.map(process_file, input_files)

