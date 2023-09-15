import math
import subprocess
import gzip
import sys
import os


# Read command-line arguments
input_filepath = sys.argv[1]
output_directory = sys.argv[2]
genome_build = sys.argv[3]

input_filename = os.path.basename(sys.argv[1])

# Create the output directory if it doesn't exist
#output_directory = os.path.join("/mnt/ebs/jackal/FILER2/FILER2-production/Homer/motif/processed_bed", genome_build)
#os.makedirs(output_directory, exist_ok=True)

output_filename = input_filename.replace(".bed.gz", ".processed.bed")
output_filepath = os.path.join(output_directory, output_filename)


# Determine the motif numbers based on input_filename
if input_filename == "FOXA1(Forkhead).bed.gz":
    motif_numbers = ["motif109", "motif110"]
    motif_length = 10
    threshold_109 = 6.337526
    threshold_110 = 6.568115
elif input_filename == "Nr5a2(NR).bed.gz":
    motif_numbers = ["motif203", "motif247"]
    motif_length = 10
    threshold_203 = 7.021126
    threshold_247 = 7.613762
elif input_filename == "OCT:OCT(POU,Homeobox).bed.gz":
    motif_numbers = ["motif259", "motif260"]
    motif_length = 15
else:
    print(f"Unsupported input file name: {input_filename}")
    sys.exit(1)

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

def calculate_motif_score(sequence, base_probabilities):
    motif_score = 0.0
    for idx, base in enumerate(sequence):
        base_index = {'a': 0, 'c': 1, 'g': 2, 't': 3}.get(base.lower())
        if base_index is not None:
            probability = base_probabilities[idx][base_index]
            score = math.log(probability / 0.25)
            motif_score += score
    return motif_score

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

def get_cellandmotif(name, parsed_dict):
    if key in parsed_dict:
        motif_id = parsed_dict[key][0]
        cell_type = parsed_dict[key][1]
        DNA_binding = parsed_dict[key][2]
        consensus_length = parsed_dict[key][3]
    else:
        motif_id = "."
        cell_type = "."
        DNA_binding = "."
        consensus_length = "."

    return motif_id, cell_type, DNA_binding, consensus_length

#  Read base probabilities from input files based on motif_numbers
base_probabilities = {}
for motif_number in motif_numbers:
    base_probabilities[motif_number] = []

    motif_file_path = f"motif_files/{motif_number}.motif"
    with open(motif_file_path, "r") as input_file:
        next(input_file)  # Skip the header line
        for line in input_file:
            probabilities = [float(p) for p in line.strip().split()]
            base_probabilities[motif_number].append(probabilities)

# Create a dictionary to store the values from parsed_subheadings.txt
parsed_dict = {}
with open("parsed_subheadings.txt", "r") as parsed_file:
    next(parsed_file)  # Skip the header line
    for line in parsed_file:
        columns = line.strip().split("\t")
        key = columns[0]
        value = (columns[0], columns[4], columns[3], columns[9])
        parsed_dict[key] = value

with gzip.open(input_filepath, "rt", encoding="utf-8") as bed_file, open(output_filepath, "w") as output_file:
#    output_file.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t")
#    output_file.write("\t".join([f"motif_score_{motif_number}" for motif_number in motif_numbers]) + "\tmotif_name\n")

    output_file.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tmotifID\tDNA_binding_domain\tcell_type\tbinding_sequence\tconsensus_sequence\t")
    output_file.write("\t".join([f"motif_score_{motif_number}" for motif_number in motif_numbers]) + "\tmotif_name\n")

    processed_intervals = set()  # To track processed intervals

    for line in bed_file:
        if line.startswith("#"):  # Skip comment lines
            continue

        columns = line.strip().split("\t")
        chrom = columns[0]
        chromStart = columns[1]
        chromEnd = columns[2]
        name = columns[3]  # Extract motif name
        strand = columns[5]
        bed_score = int(columns[4])

        # Check if the interval length matches the motif length
        interval_length = int(chromEnd) - int(chromStart) + 1
        if interval_length != motif_length:
            continue

        # Construct sequenceID using chrom, chromStart, and chromEnd
        sequenceID = f"{chrom}:{chromStart}-{chromEnd}"

        # Check if this interval has already been processed for the same motifs
        if all((chrom, chromStart, chromEnd, motif_number) in processed_intervals for motif_number in motif_numbers):
            continue

        if strand == '+':
            sequence = get_reference(genome_build, sequenceID)
        elif strand == '-':
            ref_sequence = get_reference(genome_build, sequenceID)
            sequence = get_reverse_complement(ref_sequence)

        # Calculate motif scores for the sequence using different motif files
        motif_scores = {}
        for motif_number in motif_numbers:
            motif_score = calculate_motif_score(sequence, base_probabilities[motif_number])

            # Check if the motif_number should be skipped based on the threshold
            if (motif_number == "motif203" and motif_score <= threshold_203) or (motif_number == "motif247" and motif_score <= threshold_247):
                motif_scores[motif_number] = -1  # Assign a negative value
            elif (motif_number == "motif109" and motif_score <= threshold_109) or (motif_number == "motif110" and motif_score <= threshold_110):
                motif_scores[motif_number] = -1  # Assign a negative value
            else:
                motif_scores[motif_number] = round(motif_score)

        # Determine motif names based on scores and bed score
        motif_names = [motif_number for motif_number in motif_scores if motif_scores[motif_number] == bed_score]

        if motif_names:
            for motif_name in motif_names:
                key = motif_name

                if key in parsed_dict:
                    motif_id, cell_type, DNA_binding, consensus_sequence = get_cellandmotif(motif_name, parsed_dict)
                    binding_sequence = get_reference(genome_build, sequenceID)
                    reverse_strand = get_reverse_complement(binding_sequence)

                    if strand == '+':
                        binding_sequence = get_reference(genome_build, sequenceID)
                    elif strand == '-':
                        ref_sequence = get_reference(genome_build, sequenceID)
                        binding_sequence = get_reverse_complement(ref_sequence)

                    # Write the line to the output file
                    output_line = "\t".join([
                        chrom, chromStart, chromEnd, name, str(bed_score), strand, motif_id, DNA_binding, cell_type, binding_sequence, consensus_sequence,
                        "\t".join([str(motif_scores[motif_number]) for motif_number in motif_numbers]),
                        motif_name, motif_id, cell_type, DNA_binding, consensus_sequence
                    ])
                    output_file.write(output_line + "\n")

                    # Update processed intervals
                    processed_intervals.update((chrom, chromStart, chromEnd, motif_number) for motif_number in motif_names)

print("Output written to", output_filepath)

# Compress the output file using bgzip
bgzip_output_filepath = output_filepath + ".gz"
subprocess.run(["bgzip", output_filepath])

print("Output compressed to", bgzip_output_filepath)
