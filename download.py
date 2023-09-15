import os
import requests
import gzip
import sys

def download_bed_file(url, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    response = requests.get(url)
    file_name = os.path.basename(url)
    file_path = os.path.join(output_dir, file_name)

    with open(file_path, "wb") as output_file:
        output_file.write(response.content)

    return file_path

def process_bed_file(input_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    with gzip.open(input_path, "rt", encoding="utf-8") as gz_file:
        for line in gz_file:
            if line.startswith("#"):
                continue

            columns = line.strip().split("\t")
            if len(columns) >= 4:
                fourth_column_value = columns[3]
                output_file_name = os.path.join(output_dir, f"{fourth_column_value}.bed")

                with open(output_file_name, "a") as output_file:
                    output_file.write(line)

def compress_bed_files(output_dir):
    for bed_file in os.listdir(output_dir):
        if bed_file.endswith(".bed"):
            bed_file_path = os.path.join(output_dir, bed_file)
            os.system(f"bgzip {bed_file_path}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python download_bed_files.py <genome_build>")
        sys.exit(1)

    genome_build = sys.argv[1]

    raw_files_dir = "raw_files_download"
    hg19_url = "http://homer.ucsd.edu/homer/data/motifs/homer.KnownMotifs.hg19.191020.bed.gz"
    hg38_url = "http://homer.ucsd.edu/homer/data/motifs/homer.KnownMotifs.hg38.191020.bed.gz"

    if genome_build == "hg19":
        download_url = hg19_url
    elif genome_build == "hg38":
        download_url = hg38_url
    else:
        print("Invalid genome build. Supported values: hg19, hg38")
        sys.exit(1)

    original_bed_path = download_bed_file(download_url, raw_files_dir)

    output_dir = f"split_by_motifName/{genome_build}"
    process_bed_file(original_bed_path, output_dir)
    compress_bed_files(output_dir)

    print("All BED files created and compressed.")

if __name__ == "__main__":
    main()

