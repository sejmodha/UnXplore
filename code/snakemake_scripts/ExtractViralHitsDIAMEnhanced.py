#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed  9 Oct 13:34:33 BST 2019

@author: sejmodha
"""

import argparse
import csv, os
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description='This script joins DIAMOND output file, LCA file and pileup file and, extract contigs with 100% viral hits in FASTA format')
parser.add_argument('-n','--ContigListFile', help='Input a valid contig ID list file',required=True)
parser.add_argument('-f','--FASTAFile', help='Input a valid contig file in FASTA format',required=True)
parser.add_argument('-l','--LCAFile',help='Specify a name for the LCA file generated using ExtractLCA.py',required=True)
parser.add_argument('-d','--DIAMPileupFile',help='Specify file generated using ExtractSeqStats.py (joinDIAMONDandPileup)',required=True)
parser.add_argument('-o','--outFile',help='Specify a name for the output file',required=True)
parser.add_argument('-s','--outFASTAFile',help='Specify a name for the output FASTA file with viral sequences',required=True)
args=parser.parse_args()

contiglistfile=args.ContigListFile
infasta=args.FASTAFile
inlca=args.LCAFile
diampileup=args.DIAMPileupFile
outfasta=args.outFASTAFile
outviralfile=args.outFile

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

#function to join the LCA file with pileup file and extract contigs that have only viral hits
#This function returns a list of contig IDs
def extract_viral_hits(contiglistfile,diampileupfile,lcafile,outputviralfile):
    contig_header='ContigID'
    contig_df=pd.read_csv(contiglistfile,header=None)
    contig_df.rename(columns={0:'ContigID'},inplace=True)
    #print(contig_df.head())
    diam_pileup_df=pd.read_csv(diampileupfile,sep='\t')
    lca_df=pd.read_csv(lcafile, sep='\t')
    viral_lca_df=lca_df[lca_df['ProportionOfViralHits']==100]
    contig_diam_pileup_df=diam_pileup_df[diam_pileup_df['qseqid'].isin(contig_df['ContigID'])]
    #print(contig_diam_pileup_df.head())
    viral_diam_pileup_lca_df=get_max_rows(pd.merge(contig_diam_pileup_df,viral_lca_df,left_on='qseqid',right_on='ContigID'))
    viral_diam_pileup_lca_df.to_csv(outputviralfile, index=False, header=True,sep='\t')
    viral_contigs=viral_diam_pileup_lca_df['ContigID'].tolist()
    return viral_contigs

#Function to extract given list of contig IDs from a FASTA file and save it in another FASTA file
def extract_seqs(seqid_set,fastain,fastaout):
    out_fasta=open(fastaout,'w')
    for record in SeqIO.parse(fastain, "fasta"):
        if record.id in seqid_set:
            outseq=('>'+record.id+'\n'+str(record.seq)+'\n')
            out_fasta.write(outseq)

if os.stat(contiglistfile).st_size == 0:
    os.mknod(outviralfile)
    os.mknod(outfasta)
else:
    extract_seqs(extract_viral_hits(contiglistfile,diampileup,inlca,outviralfile),infasta,outfasta)
