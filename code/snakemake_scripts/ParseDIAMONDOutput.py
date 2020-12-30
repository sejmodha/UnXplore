#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu 25 Apr 13:27:42 2019

@author: sejmodha
"""

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='This script parses BLAST tabular output')
parser.add_argument('-i','--infile', help='A valid DIAMOND (BLASTX) tabular output file',required=True)
parser.add_argument('-f','--fastafile',help='Fasta file with contigs',required=True)
parser.add_argument('-o','--outfile',help='File name for filtered list of contigs', required=True)
args=parser.parse_args()

input=args.infile
output=args.outfile
infasta=args.fastafile

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
#A function to parse DIAMOND BLAST tabular output 6 format
def extract_unmatched(blastout,fastafile,outfasta):
    blast_handle = open(blastout,'r')
    fasta_handle = SeqIO.parse(fastafile, "fasta")
    output_handle = open(output,'a')
    #Read a file line by line in a file handle
    for line in blast_handle:
        #print(line)
        fields = line.rstrip('\n').split('\t')
        seqid = fields[0]
        #If a line contains a query for which no hit is found then extract that line and
        #look for that Query ID in the fasta file and extract corresponding sequence
        if fields[1] == '*' and fields[2] == '-1':
            for record in fasta_handle:
                if seqid == record.id:
                    outseq=('>'+record.id+'\n'+str(record.seq)+'\n')
                    output_handle.write(outseq)
                    break
    blast_handle.close()
    #fasta_handle.close()
    output_handle.close()

extract_unmatched(input,infasta,output)
