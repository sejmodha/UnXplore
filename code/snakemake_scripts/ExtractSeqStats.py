#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon 9 Sep 11:03:07 BST 2019

@author: sejmodha
"""

import argparse
import pandas as pd


parser = argparse.ArgumentParser(description='This script parses pileup generated using pileup.sh (bbtools) and generates a new file that combines BLASTX (m6 - tabular) results and pipeline results')
parser.add_argument('-d','--DiamondFile', help='Input a valid DIAMOND output file in tabular (BLAST - m6)  format',required=True)
parser.add_argument('-p','--PileupFile', help='Input a valid bbmap pileup file',required=True)
parser.add_argument('-o','--outFile',help='Specify a name for the output file',required=True)
args=parser.parse_args()

indiam=args.DiamondFile
outfile=args.outFile
inpileup=args.PileupFile

def joinDIAMONDandPileup(diam_input,pileup_input,output):
    diam_header='qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qframe,staxids,stitle'.split(',')
    diam_df=pd.read_csv(diam_input,sep='\t',header=None,names=diam_header)
    #subset the dataframe to keep records with at least one hit and exclude rows that do not have a hit
    diam_df=diam_df[diam_df['sseqid'] !='*']
    #print(diam_df.shape)
    pileup_df=pd.read_csv(pileup_input,sep='\t')
    #print(pileup_df.shape)
    big_df=pd.merge(diam_df,pileup_df,left_on='qseqid',right_on='#ID')
    #print(big_df.shape)
    big_df.to_csv(output,index=False,sep='\t')

joinDIAMONDandPileup(indiam,inpileup,outfile)
