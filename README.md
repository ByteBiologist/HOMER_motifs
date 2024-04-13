# HOMER Data processing

HOMER Known Motifs 

This collection of Python scripts automates various tasks associated with processing Homer known motifs. This pipeline is designed to work with genome builds hg19/hg38 and perform motif analysis.

## Usage

Here's a overview of how to use the Homer pipeline scripts:

#### Step 1) Download and split HOMER known motifs data in BED format 
Download and split the HOMER known motifs data. Specify the genome build as an argument (hg19/hg38).
  * **python download.py hg[19/38]**

#### Step 2) Collect Motif Information
Gather information about individual motifs from the Homer database and download position frequency matrix files. This step is necessary before processing the BED files.
  * **python get_motif_info.py**

#### Step 3) Process downloaded BED Files
Process the downloaded motif files using motif information collected in step 2. Specify genome build as an argument (hg19/hg38).
  * **python generate_output_bed.py**
  
##### sub step 3a) motif score calculations
During Step 3 processing, the motif_score.py script is called upon when necessary for motif score calculations.
  * **motif_score.py**

