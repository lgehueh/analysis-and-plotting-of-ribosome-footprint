#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fixed version of plotting_reads_from_sorted_BAM_file.py
Resolves return statement bug and improves error handling
"""

import glob, re, os
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
from Bio import SeqIO, Entrez
from Bio.Seq import Seq

def gene_informn(entrez_id):
    '''
    Gets gene information from NCBI database with improved error handling
    
    Parameters
    ----------
    entrez_id : str
        The entrez id for the gene (e.g., 'NM_002300.8')

    Returns
    -------
    gene_info : dict
        Dictionary with start codon pos, stop codon pos, second in-frame 
        stop codon pos and transcript length
    '''
    Entrez.email = 'saubhiksom@iisc.ac.in'
    pattern_1 = r"\d+"
    
    try:
        with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=entrez_id) as handle:
            seq_record = SeqIO.read(handle, "gb")
    except Exception as e:
        print(f"Error fetching data for {entrez_id}: {e}")
        raise
    
    gene_info = {}
    
    for features in seq_record.features:
        if features.type == "CDS":
            transcript_seq = str(seq_record.seq)
            location = str(features.location)
            matches = re.findall(pattern_1, location)
            
            if len(matches) < 2:
                raise ValueError(f"Could not parse CDS location for {entrez_id}")
            
            CDS_location_ATG = int(matches[0])
            CDS_location_STOP = int(matches[1])
            UTR_seq = str(seq_record.seq[CDS_location_STOP:])
            ISR_protein = str(Seq.translate(seq_record.seq[CDS_location_STOP:], to_stop=True))
            
            gene_info[entrez_id] = {
                'start': CDS_location_ATG,
                'stop1': CDS_location_STOP,
                'stop2': CDS_location_STOP + len(ISR_protein) * 3 + 3,
                'length_gene': len(transcript_seq),
                'UTR_length': len(UTR_seq)
            }
            break
    
    if not gene_info:
        raise ValueError(f"No CDS feature found for {entrez_id}")
    
    return gene_info

def plotting_sorted_BAM(gene_name, entrez_id, list_of_files=None, g4_coords=None):
    '''
    Counts reads from sorted BAM files and generates plots with G4 analysis
    
    Parameters
    ----------
    gene_name : str
        Name of the gene (used for folder creation)
    entrez_id : str
        The NM_xxxxx id of the gene
    list_of_files : list, optional
        List of SRR accessions to process. If None, processes all BAM files
    g4_coords : dict, optional
        Dictionary with G4 coordinates for each entrez_id
        
    Returns
    -------
    pd.DataFrame
        DataFrame with alignment statistics and G4 region analysis
    '''
    
    # Input validation
    if not gene_name or not entrez_id:
        raise ValueError("Gene name and Entrez ID are required")
    
    if g4_coords is None:
        raise ValueError("G4 coordinates dictionary is required")
    
    if entrez_id not in g4_coords:
        raise KeyError(f"G4 coordinates not found for {entrez_id}")
    
    try:
        gene_info = gene_informn(entrez_id)
    except Exception as e:
        print(f"Error getting gene information for {entrez_id}: {e}")
        return pd.DataFrame()  # Return empty DataFrame on error
    
    # Extract gene information
    g4_start, g4_end = g4_coords[entrez_id]
    start = gene_info[entrez_id]['start']
    stop1 = gene_info[entrez_id]['stop1']
    stop2 = gene_info[entrez_id]['stop2']
    length_gene = gene_info[entrez_id]["length_gene"]
    
    # Create output directory
    location = f"./{gene_name}/"
    os.makedirs(location, exist_ok=True)
    
    # Find BAM files
    bam_files = glob.glob('sorted_bam_files//*.sorted.bam')
    
    if list_of_files is not None:
        if not isinstance(list_of_files, list):
            raise TypeError("list_of_files must be a list")
        # Filter BAM files based on provided SRR list
        bam_files = [f for f in bam_files if any(srr in f for srr in list_of_files)]
    
    if len(bam_files) == 0:
        print("No BAM files found")
        return pd.DataFrame()
    
    print(f"Processing {len(bam_files)} BAM files for {gene_name}")
    
    # Initialize data collection
    all_tables = []
    total_reads_all_srrs = 0
    a_total = np.zeros(length_gene)
    
    for filename in bam_files:
        try:
            # Extract SRR number
            srr_matches = re.findall(r'SRR\d+', filename)
            if not srr_matches:
                print(f"Could not extract SRR number from {filename}")
                continue
            srr = srr_matches[0]
            
            print(f"Processing {srr}...")
            
            # Initialize read coverage array
            a = np.zeros(length_gene)
            
            # Process BAM file
            with pysam.AlignmentFile(filename, "rb") as sorted_bamfile:
                total_reads = 0
                mapped_reads = 0
                
                for read in sorted_bamfile.fetch(until_eof=True):
                    total_reads += 1
                    
                    try:
                        if (read.reference_name == entrez_id and 
                            not read.is_unmapped):
                            
                            position = read.pos
                            length = read.qlen
                            
                            # Ensure we don't go out of bounds
                            end_pos = min(position + length, length_gene)
                            if position < length_gene and end_pos > position:
                                a[position:end_pos] += 1
                                a_total[position:end_pos] += 1
                                mapped_reads += 1
                                
                    except (ValueError, AttributeError) as e:
                        # Skip problematic reads
                        continue
                
                total_reads_all_srrs += total_reads
            
            print(f"{srr}: {mapped_reads} mapped reads out of {total_reads} total reads")
            
            # Calculate RPM (Reads Per Million)
            if total_reads > 0:
                a_rpm = (a / total_reads) * 1e6
            else:
                print(f"Warning: No reads found for {srr}")
                continue
            
            # Generate plot
            try:
                fig, ax = plt.subplots(1, 1, figsize=[16, 8], dpi=300)
                ax.plot(a_rpm)
                ax.set_xticks(ticks=[start, stop1, stop2])
                ax.set_xticklabels(labels=['S', '*', '*'])
                ax.set_ylabel('Reads per million', fontsize=18)
                ax.set_xlabel("Codon position", fontsize=18)
                fig.suptitle(f'{gene_name}_{srr}', fontsize=20)
                
                plot_filename = os.path.join(location, f'{gene_name}_{srr}.png')
                fig.savefig(plot_filename, dpi=300, bbox_inches='tight')
                plt.close(fig)  # Close figure to free memory
                
            except Exception as e:
                print(f"Error creating plot for {srr}: {e}")
            
            # Calculate G4 region statistics
            try:
                # Define regions around G4
                up_start = max(0, g4_start - 30)
                up_end = max(0, g4_start - 1)
                down_start = min(length_gene - 1, g4_end + 1)
                down_end = min(length_gene - 1, g4_end + 30)
                
                # Calculate reads per base per million (RPBM)
                if g4_end >= g4_start:
                    rd_g4 = np.sum(a[g4_start:g4_end + 1]) / (g4_end - g4_start + 1)
                else:
                    rd_g4 = 0
                
                if up_end >= up_start:
                    rd_up = np.sum(a[up_start:up_end + 1]) / (up_end - up_start + 1)
                else:
                    rd_up = 0
                
                if down_end >= down_start:
                    rd_down = np.sum(a[down_start:down_end + 1]) / (down_end - down_start + 1)
                else:
                    rd_down = 0
                
                # Create data row
                table_data = pd.DataFrame([[
                    entrez_id, gene_name, srr, length_gene, 
                    rd_up, rd_g4, rd_down, mapped_reads, total_reads
                ]], columns=[
                    "Gene ID", "Gene Name", "SRR No", "Length", 
                    "RPBM Upstream", "RPBM G4", "RPBM Downstream",
                    "Mapped Reads", "Total Reads"
                ])
                
                all_tables.append(table_data)
                
            except Exception as e:
                print(f"Error calculating G4 statistics for {srr}: {e}")
                
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            continue
    
    # Combine all results
    if all_tables:
        final_table = pd.concat(all_tables, ignore_index=True)
        print(f"\nCompleted analysis for {gene_name}")
        print(f"Processed {len(all_tables)} files successfully")
        return final_table
    else:
        print(f"No data collected for {gene_name}")
        return pd.DataFrame()
