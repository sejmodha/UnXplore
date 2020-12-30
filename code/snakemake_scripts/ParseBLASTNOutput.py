#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 26 Apr 09:44:39 BST 2019

@author: sejmodha
"""

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='This script parses BLAST tabular output')
parser.add_argument('-i','--infile', help='A valid BLAST tabular output file',required=True)
parser.add_argument('-f','--fastafile',help='Fasta file with contigs',required=True)
parser.add_argument('-o','--outfile',help='File name for filtered list of contigs', required=True)
args=parser.parse_args()

input = args.infile
output = args.outfile
infasta = args.fastafile

#Create an empty file as SNAKEMAKE rule will check if the file exists and will need this file!
output_handle = open(output,'a').close()

def is_fasta(filename):
    with open(filename,'r') as handle:
        fasta = SeqIO.parse(handle,'fasta')
        return any(fasta)

if is_fasta(infasta):
    pass
else:
    print('Input a valid FASTA file using parameter -f ')

#A function to extract Queries that do not have a BLAST hit in the BLAST output format 7
def extract_unmatched(blastout,fastafile,outfasta):
    #blast_handle=open(blastout,'r')
    fasta_handle = SeqIO.parse(fastafile, "fasta")
    output_handle = open(output, 'a')
    #Open the BLAST output file so it enable going back and forth in lines
    with open(blastout) as blast_handle:
        lines = tuple(blast_handle)
    counter=0
    for i,line in enumerate(lines):
        #Search for line that contains string indicating no hit was found
        if '# 0 hits found' in line:
            #Extract 2 lines above the no hits found line and save it in list
            query_list = list(lines[max(i-2, 0):i+1])
            #Extract the first element of the list i.e. Query ID
            query = (query_list[0].rstrip('\n').split(': ')[1])
            #print(query)
            #Open the corresponding FASTA file
            for record in fasta_handle:
                #Search for the query ID in the fasta file and extract the sequence for that ID
                if query == record.id:
                    counter=counter+1
                    outseq = ('>'+record.id+'\n'+str(record.seq)+'\n')
                    #print(outseq)
                    output_handle.write(outseq)
                    break
    print(counter)
    blast_handle.close()
    output_handle.close()

extract_unmatched(input,infasta,output)
