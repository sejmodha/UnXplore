#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 12 Jul 22:08:45 BST 2019

@author: sejmodha
"""

import argparse

parser = argparse.ArgumentParser(description='This script parses pileup generated using pileup.sh (bbtools) and generates a new pileup file for provided fasta file')
parser.add_argument('-i','--inputPileup', help='Input a valid pileup file generated using pileup.sh (BBTools)',required=True)
parser.add_argument('-o','--outputPileup',help='Output file name', required=True)
parser.add_argument('-f','--fastafile',help='FASTA file with sequences to be filtered',required=True)
parser.add_argument('-r','--runid',help='ENA/SRA run id',required=True)
parser.add_argument('-b','--bioproject',help='BioProject number',required=True)
args=parser.parse_args()


pileup_input = args.inputPileup
pileup_out = args.outputPileup
fasta_file = args.fastafile
runid = args.runid
bioproject = args.bioproject



def get_fasta_headers(fasta_file):
    fasta_headers=set()
    for line in open(fasta_file,'r'):
        if line.startswith('>'):
            line=line.rstrip('\n')
            header=line.split('>')[1]
            #print(header)
            fasta_headers.add(header)
    #convert a set to a dictionary using dict.fromkeys()
    return dict.fromkeys(fasta_headers,'')

def extract_pileup_data(input_pileup_file,header_dict,output_pileup_file,runid,bioproject):
    output=open(output_pileup_file,'w')
    output.write('#ID'+'\t'+'Avg_fold'+'\t'+'Length'+'\t'
                 +'Ref_GC'+'\t'+'Covered_percent'+'\t'
                 +'Covered_bases'+'\t'+'Plus_reads'+'\t'
                 +'Minus_reads'+'\t'+'Read_GC'+'\t'
                 +'Median_fold'+'\t'+'Std_Dev'+'\t'
                 +'Total_reads'+'\t'+'SRARunID'+'\t'
                 +'BioProject'+'\n')
    for line in open(input_pileup_file):
        line=line.rstrip('\n')
        header=runid+'_'+line.split('\t')[0]
        #print(header)
        if header in header_dict:
            plus_reads=int(line.split('\t')[6])
            minus_reads=int(line.split('\t')[7])
            total_reads=plus_reads+minus_reads
            output.write(runid+'_'+line+'\t'+str(total_reads)+'\t'+runid+'\t'+bioproject+'\n')
    output.close()


extract_pileup_data(pileup_input,get_fasta_headers(fasta_file),pileup_out,runid,bioproject)
