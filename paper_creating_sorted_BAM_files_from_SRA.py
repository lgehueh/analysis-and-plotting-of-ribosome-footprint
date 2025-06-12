#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 11:04:33 2023

@author: smelab
"""

import os
import pysam
import subprocess
from plotting_reads_from_sorted_BAM_file import plotting_sorted_BAM

# Define the SRA accession numbers and transcript filename
gene_name='sorted_bam_files'

# with open('/media/smelab/New_Volume/Debli di paper/Iwasaki_SRR_Acc_List.txt' ,'r') as f:
#     sra_accessions = f.readlines()  # Replace with your SRA accession numbers
# sra_accessions=[i.replace("\n",'') for i in sra_accessions]
#mtch2=['SRR3665720','SRR606204','SRR3623564','SRR2724718','SRR3623564','SRR1382096'] +['SRR1630813']
#ago=['SRR1605309','SRR2052945','SRR5013257','SRR2096965','SRR5345623','SRR4293695']
sra_accessions=['SRR970561', 'SRR970587', 'SRR970538','SRR970565','SRR970490','SRR970588']

transcript_filename = "transcript11.fasta"# Replace with your transcript filename

g4_coords={ "NM_002618.4": (1047,1076)} #write the G4 region span here based on trancript coordinate system
genes_to_process = [("PEX14","NM_004565.3"),
("IL4R","NM_000418.4"),
("MLL4 (KMT2B)","NM_014727.3"),
("IGFN1","NM_001164586.2")
("TRAF5","NM_001033910.3")
("ZFYVE19","NM_001077268.2")
("ZNF687","NM_020832.3")
("ZNF830","NM_052857.4")
("SUZ12","NM_015355.4")
("FN3KRP","NM_024619.4")
("FMR1","NM_002024.6")]

# Build the bowtie2 index for the transcript
bowtie2_build_command = ["/home/smelab/bowtie2/bowtie2-build", transcript_filename, "transcript"]
subprocess.check_call(bowtie2_build_command)
print("Bowtie2 index buiding done")
# Download SRA and Align the reads to the transcript
aligned_reads = []
total_reads=0
for sra_accession in sra_accessions:
    fastq_filename ='/home/smelab/'+ f"{sra_accession}.fastq"
    if os.path.isfile(fastq_filename)==False:
        print('starting download')
        fastq_dump_command = ["/home/smelab/sratoolkit.3.0.5-ubuntu64/bin/fastq-dump", sra_accession,'--outdir','/home/smelab/']
        subprocess.check_call(fastq_dump_command)
        print('download fastq done')
    '''
    print('trimming started')
    trimmed_fastq_filename = f"{sra_accession}_trimmed.fastq"
    fastp_command = [
        "fastp",  # Assuming fastp is now accessible from the Conda environment
        "--in1", fastq_filename,
        "--out1", trimmed_fastq_filename,
        "--trim_poly_x",  # Maximum number of 'N' bases allowed in a read
        "--trim_front1", "5",  # Trim 5 bases from the 5' end (change the value as needed)
        '--thread', '4',      # worker thread number, default is 3 (int [=3])
        # Add more parameters as needed (e.g., quality filtering, length trimming)
    ]
    subprocess.check_call(fastp_command)
    print('trimming done')
    os.remove(fastq_filename)
    '''

    sam_filename = f"{sra_accession}_aligned_reads.sam"
    bowtie2_command = ["/home/smelab/bowtie2/bowtie2", "-x", 
                       "transcript", "-U", 
                       fastq_filename, "-S", 
                       sam_filename, "-p", "1", "--very-sensitive-local"]
    try:
        subprocess.run(bowtie2_command, check=True)
        # os.remove(fastq_filename)
    except subprocess.CalledProcessError as e:
        print(f"Error running Bowtie 2: {e}")
        # os.remove(fastq_filename)
        continue
    
    # Convert the SAM file to a BAM file
    bam_filename = f"{sra_accession}_aligned_reads.bam"
    samfile = pysam.AlignmentFile(sam_filename, "r")
    bamfile = pysam.AlignmentFile(bam_filename, "wb", template=samfile)
    for read in samfile:
        bamfile.write(read)
    samfile.close()
    bamfile.close()
    
    # Sort the BAM file
    sorted_bam_filename = gene_name+'/'+f"{sra_accession}_aligned_reads.sorted.bam"
    pysam.sort("-o", sorted_bam_filename, bam_filename)
    os.remove(bam_filename)

    # Index the sorted BAM file
    pysam.index(sorted_bam_filename)

    # Delete intermediate files
    os.remove(sam_filename)
      
    # Get alignment information
    aligned_reads.append((sra_accession, sorted_bam_filename))

    # Print the number of ['SRR1257257','SRR1630813']sra_accessions=mapped reads for each SRA
    sorted_bamfile = pysam.AlignmentFile(sorted_bam_filename, "rb")
    num_mapped_reads = 0
    for read in sorted_bamfile:
        if not read.is_unmapped:
            num_mapped_reads += 1
    sorted_bamfile.close()
    print(f"Number of mapped reads for {sra_accession}: {num_mapped_reads}")
    total_reads+=num_mapped_reads
    


# Print the total number of mapped reads across all SRAs
print(f"Total number of mapped reads: {total_reads}")

from plotting_reads_from_sorted_BAM_file import plotting_sorted_BAM

all_gene_data = []

for gene_name, entrez_id in genes_to_process:
    print(f"\nProcessing {gene_name} ({entrez_id})...")
    if entrez_id not in g4_coords:
        print(f"G4 coordinates not defined for {entrez_id}. Skipping.")
        continue
    df = plotting_sorted_BAM(
        gene_name,
        entrez_id,
        list_of_files=sra_accessions,
        g4_coords=g4_coords  # ← pass G4 positions
    )
    all_gene_data.append(df)

# Merge and save all results
if all_gene_data:
    final_df = pd.concat(all_gene_data, ignore_index=True)
    final_df.to_excel("G4_RPBM_All_Genes.xlsx", index=False)
    print("Saved: G4_RPBM_All_Genes.xlsx")
else:
    print("No gene data was processed.")


# Count the number of index files
index_files = [filename for filename in os.listdir(".") if filename.startswith("transcript.")]
num_index_files = len(index_files)
print(f"Number of index files: {num_index_files}")

# Clean up transcript index files
for index_file in index_files:
    os.remove(index_file)

# AQP=plotting_sorted_BAM('AQP4','NM_001364286.1',list_of_files=sra_accessions)
# VEGF=plotting_sorted_BAM('VEGFA', 'NM_001025366.3',list_of_files=sra_accessions)
# MAPK10=plotting_sorted_BAM('MAPK10', 'NM_138982.4',list_of_files=sra_accessions)
# plotting_sorted_BAM('LDHB','NM_002300.8')
# d1=plotting_sorted_BAM('MTCH2', 'NM_019758.3',list_of_files=sra_accessions)
