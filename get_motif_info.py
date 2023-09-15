import re
import os
import requests
from bs4 import BeautifulSoup

def get_first_subheading(url):
    try:
        # Send a GET request to the URL to fetch the HTML content
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for unsuccessful requests

        # Parse the HTML content using BeautifulSoup
        soup = BeautifulSoup(response.content, 'html.parser')

        # Find the first h2 element (sub-heading) in the HTML page
        h2_element = soup.find('h2')

        if h2_element:
            # Store the content of the first sub-heading
            subheading_content = h2_element.text.strip()
            return subheading_content
        else:
            print(f"No sub-heading (h2 element) found on the webpage: {url}")
            return None
    except Exception as e:
        print(f"Error occurred while processing URL {url}: {e}")
        return None

def get_motif_sequence_and_length(motif_name):
    base_url = f'http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults/{motif_name}.motif'
    response = requests.get(base_url)

    text_content = response.text
    match = re.search(r">(\S+)", text_content)

    if match:
        sequence = match.group(1)
        motif_length = len(sequence)
    else:
        sequence = "N/A"
        motif_length = 0

    return sequence, motif_length

def parse_subheading(line):
    parts = line.strip().split()
    if len(parts) >= 5:

        # Get the motif ID from the fifth column
        motif_name = "motif" + parts[4].replace("(", "").replace(")", "").lower() 

        # Third column (AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer)
        third_column = parts[2]

        # Split third_column by '/' into 3 parts
        split_columns = third_column.split('/')

        # adding name column
        name = split_columns[0]

        # Parse the first split column into TF-name and DNA-binding-domain
        missing_DNA_binding = "ZNF652" # Dealing with missing DNA binding domain info
        misspelled_tf_name = "LHX9" # Dealing with misspelled TF name

        if missing_DNA_binding in split_columns[0]:
            tf_name = split_columns[0]
            dna_domain = '.'
        elif misspelled_tf_name in split_columns[0]:
            tf_name, dna_domain = split_columns[0].split('(')
            tf_name = "LXH9"
            name = "LXH9(Homeobox)"
        else:
            tf_name, dna_domain = split_columns[0].split('(')
            dna_domain = dna_domain.replace(')', '').strip()
            tf_name = tf_name.strip()

        # Parse the second split column ThioMac-PU.1-ChIP-Seq(GSE21512) into ThioMac \t PU.1 \t ChIP-Seq \t GSE21512
        cell_type, immunoprecipitated_protein, assay, seq = "", "", "", "" # Initialize variables with empty strings

        tmp = split_columns[1]

        # Handling Expressions
        if 'Expression' in tmp:
            cell_type, immunoprecipitated_protein, assay = tmp.split('-', 2)
            assay = re.sub(r'\([^)]*\)', '', assay)  # Remove anything within parentheses, ThioMac-LPS-Expression(GSE23622)
            # Some expression motifs have origin info, ThioMac-LPS-Expression(GSE23622)
            if '(' in tmp:
                origin = tmp.split('(')[1].rstrip(')')
            else:
                origin = "."

        # Handling promoters 
        elif 'Promoter' in tmp:
            cell_type = "."
            immunoprecipitated_protein = "."
            assay = tmp
            origin = "."
        else:

            # Some motifs do not have origin info K562-JunD-ChIP-Seq()
            try:
                origin = tmp.split('(')[1].rstrip(')')
            except IndexError:
                origin = "."

            if misspelled_tf_name in tmp: # Replacing misspelled TF name with the correct name
                tmp = "Hct116-LXH9.V5-ChIP-Seq"

            # Extracting info ZebrafishEmbryos-Cdx4.Myc-ChIP-Seq(GSE48254)
            if len(tmp.split('-')) >= 4:
                cell_type, immunoprecipitated_protein, assay, seq = tmp.split('-', 3)
                seq = re.sub(r'\([^)]*\)', '', seq)  # Remove anything within parentheses, mES-cMyc-ChIP-Seq(GSE11431)
                assay = assay + '-' + seq

            # Extracting info HUDEP2-KLF1-CutnRun(GSE136251)
            if len(tmp.split('-')) == 3:
               split_parts = tmp.split('-', 2)
               # Handling K562-mStart-Seq (has 3 parts but no immunoprecipitated_protein)
               if "mStart" in split_parts[1]:
                   cell_type, assay, seq = split_parts
                   assay = assay + '-' + seq
                   immunoprecipitated_protein = "."
               else:
                   cell_type, immunoprecipitated_protein, assay = tmp.split('-', 2)
                   assay = re.sub(r'\([^)]*\)', '', assay)  # Remove anything within parentheses, mES-cMyc-ChIP-Seq(GSE11431)

        resource = split_columns[-1].strip()

        consensus_sequence, motif_length = get_motif_sequence_and_length(motif_name)

        # Create a list containing the parsed values in desired order
        output_list = [ motif_name, name, tf_name, dna_domain, cell_type, immunoprecipitated_protein, assay, origin, resource, consensus_sequence, str(motif_length)]

        return "\t".join(output_list)
    else:
        return None

def download_motif_file(url, save_path, file_name):
    try:
        # Send a GET request to the URL to get the file contents
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for unsuccessful requests

        # Save the contents to a local file with the specified file name
        file_path = os.path.join(save_path, file_name)
        with open(file_path, 'wb') as f:
            f.write(response.content)

        print(f"File {file_name} downloaded and saved successfully.")
        return True
    except Exception as e:
        print(f"Error occurred: {e}")
        return False


if __name__ == "__main__":
    base_url = "http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults/motif"
    num_motifs = 436

    with open("subheadings.txt", "w") as file:
        for i in range(1, num_motifs + 1):
            url = f"{base_url}{i}.info.html"
            motif_subheading = get_first_subheading(url)
            if motif_subheading:
#                file.write(f"Sub-heading for {url}:\n")
                file.write(motif_subheading + "\n")

    with open("parsed_subheadings.txt", "w") as output_file:
        header = "motif_name\tname\ttf_name\tDNA_binding_domain\tcell_type\timmunoprecipitated_protein\tassay\torigin\tresource\tconsensus_sequence\tmotif_length\n"
        output_file.write(header)

        with open("subheadings.txt", "r") as input_file:
            for line in input_file:
                parsed_line = parse_subheading(line)
                if parsed_line:
                    output_file.write(parsed_line + "\n")

    # Download motif matrix files
    matrix_base_url = "http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults/motif"
    save_matrix_directory = "motif_files"

    if not os.path.exists(save_matrix_directory):
        os.makedirs(save_matrix_directory)

    num_matrix_files = num_motifs  # Same number of matrix files as motif files

    for i in range(1, num_matrix_files + 1):
        matrix_file_name = f"motif{i}.motif"
        matrix_file_url = f"{matrix_base_url}{i}.motif"
        download_motif_file(matrix_file_url, save_matrix_directory, matrix_file_name)

