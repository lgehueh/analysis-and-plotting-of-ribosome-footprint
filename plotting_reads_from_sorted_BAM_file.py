#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 17:54:08 2023

@author: smelab
"""

import glob, re,os
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
# from Bio.Blast import NCBIXML



def gene_informn(entrez_id):
    '''
    This code gets information for a given entrez id from NCBI database. Stable
    internet connection is needed for this code

    Parameters
    ----------
    entrez_id : String
        The entrez id for the gene.

    Returns
    -------
    gene_info : Dictionary
        Return Start codon pos, stop codon pos, second in-frame stop codon pos
        and transcript length for the input entrez id. Uncomment the lines for
        other details if needed

    '''
    Entrez.email = 'saubhiksom@iisc.ac.in'
    pattern_1="\d+"
    with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=entrez_id) as handle:
        seq_record = SeqIO.read(handle, "gb") # using "gb" as an alias for "genbank"
    gene_info={}               
    for features in seq_record.features:
        if features.type == "CDS":
            transcript_seq=str(seq_record.seq)
            location=str(features.location)
            CDS_location_ATG= int(re.findall(pattern_1, location)[0])
            CDS_location_STOP=int(re.findall(pattern_1, location)[1])
            # gene_name=seq_record.description
            # CDS_seq=str(seq_record.seq[CDS_location_ATG:CDS_location_STOP])
            UTR_seq=str(seq_record.seq[CDS_location_STOP:])
            ISR_protein=str(Seq.translate(seq_record.seq[CDS_location_STOP:], to_stop=True))
            # ISR_seq=str(seq_record.seq[CDS_location_STOP:CDS_location_STOP+len(ISR_protein)*3+3])
            # Rem_UTR_seq=str(seq_record.seq[CDS_location_STOP+len(ISR_protein)*3+3:])
            # assert len(ISR_seq)+len(Rem_UTR_seq)==len(UTR_seq), "position error"
            gene_info[entrez_id]={'start':CDS_location_ATG,
                                                         "stop1":CDS_location_STOP,
                                                         "stop2":CDS_location_STOP+len(ISR_protein)*3+3,
                                                         "length_gene":len(transcript_seq),
                                                         "UTR_length":len(UTR_seq)}
            break
    return gene_info
def plotting_sorted_BAM(gene_name,entrez_id,list_of_files=None):
    '''
    This function counts the reads from a sorted BAM file and returns a plot.

    Parameters
    ----------
    gene_name : Name of the gene. A folder will be created in the working 
                directory with same name
    entrez_id : The NM_xxxxx id of the gene
    
    list_of_files: Give list of SRRs if you want. Those fies must be in sorted
                    BAM file folder

    Returns:
        Plots the numbers of reads per million for the input transcript
        returns a pandas dataframe with alignment info in CDS and UTR
    -------
    '''
    gene_info=gene_informn(entrez_id)
    start = gene_info[entrez_id]['start']
    stop1 = gene_info[entrez_id]['stop1']
    stop2 = gene_info[entrez_id]['stop2']
    length_gene=gene_info[entrez_id]["length_gene"]
    utr_length=gene_info[entrez_id]["UTR_length"]

    location = "./%s/" %gene_name
    if os.path.isdir(location) == False:
        os.makedirs(location)
    bam_files=glob.glob('sorted_bam_files//'+"*.sorted.bam")
    if list_of_files != None:
        if type(list_of_files)!= list:
            print("List input needed")
            return
        else:bam_files=[i for i in bam_files if (re.findall('SRR\d*', i)[0]) in list_of_files]
    if len(bam_files)==0:
        print("no bam files found")
        return
    total_reads_all_srrs=0
    a_total=np.zeros(length_gene)
    # table_data = pd.DataFrame(columns=["gene id", "gene name", "SRR no", "transcript/n length", "RD CDS", "RD 3'-ISR"])
    for filename in bam_files:
        srr=re.findall('SRR\d*', filename)[0]
        a=np.zeros(length_gene)
        sorted_bamfile = pysam.AlignmentFile(filename,"rb")
        total_reads=0
        for read in sorted_bamfile.fetch(until_eof=True):
            total_reads+=1
            total_reads_all_srrs+=1
            try:
                if read.reference_name == entrez_id:# and read.qlen >= 27: # and 3<=len(read.cigarstring)<=4: #and read.qlen < 30:
                    # print(read.cigarstring)
                    position = read.pos
                    length = read.qlen
                    a[position - 1 : position - 1 + length] += 1
                    a_total[position - 1 : position - 1 + length] += 1
                    # print(read.reference_name, read.cigarstring, position, length)
            except ValueError:
                continue

        sorted_bamfile.close()        
        
        a_rpm = (a / total_reads) * 1e6 #calculating RPM
        fig,ax=plt.subplots(1,1,figsize=[16,8],dpi=600)
        ax.plot(a_rpm)
        ax.set_xticks(ticks=[start,stop1, stop2])
        ax.set_xticklabels(labels=['S','*','*',])
        # plt.xlim(0,stop2+1000)
        ax.set_ylabel('Reads per million', fontsize=18)
        ax.set_xlabel("Codon position", fontsize=18)
        fig.suptitle(gene_name+'_'+re.findall('SRR\d*', filename)[0], fontsize=20 )
        fig.savefig(location+gene_name+'_'+re.findall('SRR\d*', filename)[0]+'.png', dpi=300)
        plt.show()
        RD_CDS=np.sum(a[start+12:stop1+1])/len(a[start+12:stop1+1])
        RD_ISR=np.sum(a[stop1+12:stop2+1])/len(a[stop1+12:stop2+1])
        RD_rem_utr=np.sum(a[stop2+12:stop2+212])/len(a[stop2+12:stop2+212])
        percent_RT=RD_ISR/RD_CDS*100
        table_data=np.array([entrez_id,gene_name ,srr,length_gene,utr_length,RD_CDS ,RD_ISR,RD_rem_utr,percent_RT])
        table_data=pd.DataFrame(table_data.reshape(1,len(table_data)),columns=["gene id","gene name","SRR no","transcript\nlength","3'-UTR\nlength",'RD CDS',"RD ISR","RD rem 3'-UTR","percent RT"])
    '''        
    a_total_rpm = (a_total / total_reads_all_srrs ) * 1e6 #calculating RPM
    fig,ax=plt.subplots(1,1,figsize=[16,8],)
    ax.plot(a_total_rpm)
    ax.set_xticks(ticks=[start,stop1, stop2])
    ax.set_xticklabels(labels=['S','*','*',])
    # plt.xlim(0,stop2+1000)
    ax.set_ylabel('Combined reads per million', fontsize=18)
    ax.set_xlabel("Codon position", fontsize=18)
    fig.suptitle("Combined reads", fontsize=20 )
    # fig.savefig(location+re.findall('SRR\d*', filename)[0]+'.png', dpi=300)
    plt.show()
    '''
    return table_data

# with open("file.txt","a+") as file1:
# 	for ele in a:
# 		file1.write(str(ele)+"\n")
#break
# d1=plotting_sorted_BAM('LDHB','NM_002300.8',['SRR1630833'])
# d11=plotting_sorted_BAM('LDHB','NM_002300.8',['SRR1630831'])
# d2=plotting_sorted_BAM('MDH1', 'NM_005917.4',['SRR1630831'])
# d22=plotting_sorted_BAM('MDH1', 'NM_005917.4',['SRR1630833'])
# # d3=plotting_sorted_BAM('AGO1', 'NM_012199.5',['SRR1257257'])
# # d4=plotting_sorted_BAM('MTCH2','NM_014342.4', ['SRR1630830'])
# excel_data=pd.concat([d1,d11,d2,d22],axis=0,ignore_index=True)
# excel_data.to_excel('ldh_mdh_excel_data.xlsx',index=False)