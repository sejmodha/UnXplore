#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed  9 Oct 13:34:33 BST 2019

@author: sejmodha
"""

import argparse
import csv
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description='This script joins LCA file and pileup file and extract contigs with 100% viral hits in FASTA format')
parser.add_argument('-p','--PileupFile', help='Input a valid Pileup file (from bbmap)',required=True)
parser.add_argument('-f','--FASTAFile', help='Input a valid contig file in FASTA format',required=True)
parser.add_argument('-l','--LCAFile',help='Specify a name for the LCA file generated using ExtractLCA.py',required=True)
parser.add_argument('-o','--outFile',help='Specify a name for the output file that would combine LCA and Pileup info',required=True)
parser.add_argument('-s','--outFASTAFile',help='Specify a name for the output FASTA file with viral sequences',required=True)
args=parser.parse_args()

pileup=args.PileupFile
infasta=args.FASTAFile
inlca=args.LCAFile
outfasta=args.outFASTAFile
outfile=args.outFile

#Function to check whether input is fasta
def is_fasta(filename):
    with open(filename,'r') as handle:
        fasta = SeqIO.parse(handle,'fasta')
        return any(fasta)

if is_fasta(infasta):
    pass
else:
    print('Input a valid FASTA file using parameter -f \n')
#function to join the LCA file with pileup file and extract contigs that have only viral hits
#This function returns a list of contig IDs
def extract_viral_hits(lcafile,pileupfile,lca_pileup_out):
    lcadf=pd.read_csv(lcafile,sep='\t')
    viral_lca_df=lcadf[lcadf['ProportionOfViralHits']==100]
    pileup_df=pd.read_csv(pileupfile,sep='\t')
    viral_lca_pileup=pd.merge(viral_lca_df,pileup_df,left_on='ContigID',right_on='#ID')
    viral_lca_pileup.to_csv(lca_pileup_out,index=False)
    viral_contigs=viral_lca_pileup['ContigID'].tolist()
    return viral_contigs
    
#Function to extract given list of contig IDs from a FASTA file and save it in another FASTA file
def extract_seqs(seqid_set,fastain,fastaout):
    out_fasta=open(fastaout,'w')
    for record in SeqIO.parse(fastain, "fasta"):
        if record.id in seqid_set:
            outseq=('>'+record.id+'\n'+str(record.seq)+'\n')
            out_fasta.write(outseq)


extract_seqs(extract_viral_hits(inlca,pileup,outfile),infasta,outfasta)
