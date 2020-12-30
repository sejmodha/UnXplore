#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon 12 Aug 11:44:29 BST 2019

@author: sejmodha
"""

import argparse
import csv
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description='This script parses pileup generated using pileup.sh (bbtools) and generates a new pileup file for provided fasta file')
parser.add_argument('-i','--inFile', help='Input a valid DIAMOND output file in tabular (BLAST - m6)  format',required=True)
parser.add_argument('-f','--fastaFile', help='Input a valid contig file in FASTA format',required=True)
parser.add_argument('-o','--outFastaFile',help='Specify a name for the output FASTA file',required=True)
parser.add_argument('-s','--outDiamondFile',help='Specify a name for the output DIAMOND file',required=True)
parser.add_argument('-n','--outContigNames',help='Specify a namae for the output contig list file',required=True)
args=parser.parse_args()

indiam=args.inFile
outdiam=args.outDiamondFile
infasta=args.fastaFile
outfasta=args.outFastaFile
outnames=args.outContigNames

#Function to check whether input is fasta
def is_fasta(filename):
    with open(filename,'r') as handle:
        fasta = SeqIO.parse(handle,'fasta')
        return any(fasta)

if is_fasta(infasta):
    pass
else:
    print('Input a valid FASTA file using parameter -f \n')

#Function to group rows by query seqid and get the max of the percent identity for each query
def get_max_rows(df):
    B_maxes = df.groupby('qseqid').pident.transform(max)
    return df[df.pident == B_maxes]

def extract_hits(inputfile,outputfile,contiglistfile,percid):
    '''Function to extract DIAMOND hits with <= specified percent identity threshold.
    Returns a unique set of query IDs that fulfil the criteria.
    It takes 3 arguments: input file (DIAMOND output in tabular format) output file (subset of hits that fall within the criteria) and percent identity threshold'''
    diam_header='qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qframe,staxids,stitle'.split(',')
    df=pd.read_csv(inputfile,sep='\t',header=None,names=diam_header)
    #Exclude rows where pident could be -1 as those rows represent contigs without a BLASTX hit
    df=df[(df['pident'] > 0)]
    subset_df=get_max_rows(df)
    seqid_set=set(subset_df[(subset_df['pident'] > percid)]['qseqid'].tolist())
    print(len(seqid_set))
    df[df['qseqid'].isin(seqid_set)].to_csv(outputfile,index=False,sep='\t',header=None)
    contigout=open(contiglistfile,'w')
    contigout.write('\n'.join(seqid_set))
    contigout.close()
    print(sorted(seqid_set))
    return seqid_set

def extract_seqs(seqid_set,fastain,fastaout):
    out_fasta=open(fastaout,'w')
    for record in SeqIO.parse(fastain, "fasta"):
        if record.id in seqid_set:
            outseq=('>'+record.id+'\n'+str(record.seq)+'\n')
            out_fasta.write(outseq)

extract_seqs(extract_hits(indiam,outdiam,outnames,percid=80),infasta,outfasta)
