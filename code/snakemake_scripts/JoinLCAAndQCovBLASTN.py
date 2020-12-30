#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon 14 Oct 16:46:50 BST 2019

@author: sejmodha
"""

import pandas as pd
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='This script joins the output generated using ExtractLCABLASTM6 with BLAST tabular output file generates a consolidated output file')
parser.add_argument('-l','--infile', help='Input a valid file generated using ExtractLCABLASTM6.py script',required=True)
parser.add_argument('-b','--BLASTNfile', help='Input a valid BLAST input file with following columns \nqseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qframe,staxids,stitle,qcovs ',required=True)
parser.add_argument('-p','--pileupfile',help='Input pileup file generated using bbtools pileup.sh',required=True)
parser.add_argument('-f','--fastafile',help='Input FASTA file with contigs that were used to generate supplied files',required=False)
parser.add_argument('-o','--outfile',help='Output file name', required=True)


args=parser.parse_args()

inputfile = args.infile
blastfile=args.BLASTNfile
outputfile = args.outfile
pileup=args.pileupfile
fasta=args.fastafile
#Function to group rows by query seqid and get the max of the percent identity for each query
def get_max_rows(df,colname):
    #print(colname)
    B_maxes = df.groupby('qseqid')[colname].transform(max)
    return df[df[colname] == B_maxes]

def extract_fasta_headers(fastafile):
    headers_list=[]
    for record in SeqIO.parse(fastafile,'fasta'):
        headers_list.append(record.id)
    return headers_list

blast_header='qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qframe,staxids,stitle,qcovs'
blastn_df=pd.read_csv(blastfile,names=blast_header.split(','),sep='\t')
#lca_header='ContigID,LCA_TaxonID,LCA_TaxonName,NumberOfUniqueTaxIDs,NumberOfViralHits,ProportionOfViralHits,LCA_Family,LCA_Superkingdom'
lca_df=pd.read_csv(inputfile,sep='\t')

blastn_df=get_max_rows(blastn_df,'qcovs')
blastn_df=get_max_rows(blastn_df,'sseqid')
blastn_df=blastn_df.drop_duplicates(['qseqid'])

pileup_df=pd.read_csv(pileup,sep='\t')

lca_blast_df=pd.merge(lca_df,blastn_df[['qseqid','qcovs','pident']],left_on='ContigID',right_on='qseqid')


if fasta:
    no_blastn_hit_list=list(set(extract_fasta_headers(fasta))-set(lca_blast_df.ContigID))
    #print(no_blastn_hit_list)
    for each_contig in no_blastn_hit_list:
        lca_blast_df=lca_blast_df.append({'ContigID':each_contig,'LCA_Superkingdom':'NoHit'},ignore_index=True)
    big_df=pd.merge(lca_blast_df,pileup_df,left_on='ContigID',right_on='#ID')
else:
    big_df=pd.merge(lca_blast_df,pileup_df,left_on='ContigID',right_on='#ID')

big_df.to_csv(outputfile,index=False)
