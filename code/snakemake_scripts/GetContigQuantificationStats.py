#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 12 Jul 22:08:45 BST 2019

@author: sejmodha
"""

import argparse

parser = argparse.ArgumentParser(description='This script parses pileup generated using pileup.sh (bbtools) and generates a new pileup file for provided fasta file')
parser.add_argument('-c','--contigsFile', help='Input a valid contig file in FASTA format',required=True)
parser.add_argument('-l','--longContigFile',help='Input a valid long contigs (>300nt) file in FASTA format', required=True)
parser.add_argument('-u','--unknownContigFile',help='Input a valid dark contigs file in FASTA format',required=True)
parser.add_argument('-g','--greyContigFile',help='Input a valid grey contigs file in FASTA format',required=True)
parser.add_argument('-n','--blastnOutputFile',help='Input a BLAST output file - outfmt 7 with comment lines',required=True)
parser.add_argument('-d','--diamondOutputFile',help='Input a DIAMOND output file - including records without hits',required=True)
parser.add_argument('-r','--runid',help='ENA/SRA run id',required=True)
parser.add_argument('-b','--bioproject',help='BioProject number',required=True)
parser.add_argument('-o','--outputFile',help='Specify a name for the output file',required=True)
args=parser.parse_args()

rawcontigfile=args.contigsFile
long_contigs_file=args.longContigFile
grey_contigs_file=args.greyContigFile
dark_contigs_file=args.unknownContigFile
diam_output_file=args.diamondOutputFile
blastn_output_file=args.blastnOutputFile
runid=args.runid
bioproject=args.bioproject
output=args.outputFile

#function to count number of records in a fasta file
def count_FASTA_headers(fastafile):
    contigs=0
    fh = open(fastafile,'r')
    for line in fh:
        if line.startswith(">"):
            contigs += 1
        #print(long_contigs)
    fh.close()
    return contigs

#This function parses a file in BLASTX tabular - M8 output format generated using --unal
def parse_DIAMOND_M8_output(diamond_output_file):
    diam_hits=[]
    no_diam_hits=0
    diamond_handle=open(diamond_output_file,'r')
    #Read a file line by line in a file handle
    for line in diamond_handle:
        fields=line.rstrip('\n').split('\t')
        #print(fields)
        seqid=fields[0]
        #If a line contains a query for which no hit is found then extract that line
        if fields[1] == '*' and fields[2] == '-1':
            no_diam_hits=no_diam_hits+1
        else:
            diam_hits.append(seqid)
    diamond_handle.close()
    #print(str(len(list(set(diam_hits)))))
    return str(len(list(set(diam_hits))))+','+str(no_diam_hits)

#Function to count number of records with and without a BLAST hit
#This function parses a file in BLAST tabular M7 with comments output format
def parse_BLAST_M7_output(blast_m7_output_file):
    blastn_hits=[]
    blastn_handle=open(blast_m7_output_file,'r')
    for line in blastn_handle:
        if line.startswith('#'):
            pass
        else:
            fields=line.rstrip('\n').split('\t')
            seqid=fields[0]
            blastn_hits.append(seqid)
    blastn_handle.close()
    return str(len(list(set(blastn_hits))))

raw_contigs=count_FASTA_headers(rawcontigfile)
long_contigs=count_FASTA_headers(long_contigs_file)
grey_contigs=count_FASTA_headers(grey_contigs_file)
contigs_unknown_origin=count_FASTA_headers(dark_contigs_file)
diam_stats=parse_DIAMOND_M8_output(diam_output_file)
contigs_with_diam_hits=diam_stats.split(',')[0]
contigs_without_diam_hits=diam_stats.split(',')[1]
contigs_with_blastn_hits=parse_BLAST_M7_output(blastn_output_file)
if float(contigs_unknown_origin) == 0 or float(long_contigs) == 0:
    unknown_contig_prop=0
else:
    unknown_contig_prop=round((float(contigs_unknown_origin)*100)/float(long_contigs),2)

outfile=open(output,'w')
outfile.write('BioProject'+','
              +'SRARunID'+','
              +'TotalContigs'+','
              +'LongContigs(>=300)'+','
              +'ContigsWithDIAMONDHit'+','+'ContigsWithoutDIAMONDHit'+','
              +'ContigsWithBLASTNHit'+','
              +'ContigsWithPartialHit(Grey)'+','
              +'ContigsWithUnknownOrigin'+','
              +'UnknownContigsProportion'+'\n')

outfile.write(bioproject+','
          +runid+','
          +str(raw_contigs)+','
          +str(long_contigs)+','
          +str(contigs_with_diam_hits)+','+str(contigs_without_diam_hits)+','
          +str(contigs_with_blastn_hits)+','
          +str(grey_contigs)+','
          +str(contigs_unknown_origin)+','
          +str(unknown_contig_prop)+'\n')
outfile.close()
