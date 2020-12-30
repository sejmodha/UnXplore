#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 16:29:05 2019

@author: sejmodha
"""

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='This script filters FASTA files of contigs')
parser.add_argument('-i','--infile', help='Input a valid fasta file',required=True)
parser.add_argument('-o','--outfile',help='Output a file name for filtered fasta file', required=True)
parser.add_argument('-l','--length',help='Specify a length to filter sequences', required=True)
parser.add_argument('-p','--prefix',help='Specify a prefix to rename sequences', required=False)
args=parser.parse_args()

infasta = args.infile
outfasta = args.outfile
length = args.length
pref = args.prefix

#Create an empty file as SNAKEMAKE rule will check if the file exists and will need this file!
output_handle = open(outfasta,'a').close()

try:
    val=int(length)
except ValueError:
    print('Enter a valid integer for length -l parameter')

def is_fasta(filename):
    with open(filename,'r') as handle:
        fasta = SeqIO.parse(handle,'fasta')
        return any(fasta)

if is_fasta(infasta):
    pass
else:
    print('Input a valid fasta file using parameter -i ')

if pref is None:
    pref=''

def filter_fasta_by_len(input,output,length,prefix=''):
    if is_fasta(infasta):
        long_seqs=[]
        input_handle = open(input,'r')
        output_handle = open(output,'w')
        for record in SeqIO.parse(input_handle, "fasta") :
            if len(record.seq) >= int(length) :
                #print('>'+str(prefix)+'_'+record.id+'\n'+str(record.seq))
                outseq='>'+str(prefix)+'_'+record.id+'\n'+str(record.seq)+'\n'
                output_handle.write(outseq)

filter_fasta_by_len(infasta,outfasta,length,pref)
