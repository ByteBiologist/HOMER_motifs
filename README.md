# HOMER Data processing

This repository contains a set of Python scripts to automate various tasks related to processing Homer known motifs. The pipeline is designed to work with different genome builds and perform motif analysis.

Usage
Here's a general overview of how to use the Homer pipeline scripts:

1. Download Raw BED Files: Use **download.py** to download the raw BED files. Provide the desired genome build as an argument [hg19, hg38]. This step sets up the initial dataset.

2. Collect Motif Information (One-Time Step): Run **get_motif_info.py** to gather information about individual motifs from the Homer database. This step is necessary before processing the BED files.

3. Process BED Files: Execute **generate_output_bed.py** to process the downloaded BED files using the collected motif information. Provide the genome build as an argument [hg19, hg38].

4. Advanced Processing: Use **motif_score.py** for additional processing steps, such as motif score calculations, on specific files. Ensure that the necessary motif information has been collected in step 2.

