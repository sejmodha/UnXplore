#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed 14 Aug 10:51:07 BST 2019

@author: sejmodha
"""

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='This script calculates length of given fasta file in bases')
parser.add_argument('-c','--contigsFile', help='Input a valid contig file in FASTA format',required=True)
parser.add_argument('-l','--longContigFile',help='Input a valid long contigs (>300nt) file in FASTA format', required=True)
parser.add_argument('-g','--PartiallyKnownContigFile',help='Input a valid grey contigs (>300nt) file in FASTA format', required=True)
parser.add_argument('-u','--unknownContigFile',help='Input a valid dark contigs file in FASTA format',required=True)
parser.add_argument('-r','--runid',help='ENA/SRA run id',required=True)
parser.add_argument('-b','--bioproject',help='BioProject number',required=True)
parser.add_argument('-o','--outputFile',help='Specify a name for the output file',required=True)
args=parser.parse_args()


rawcontigfile=args.contigsFile
long_contigs_file=args.longContigFile
grey_contigs_file=args.PartiallyKnownContigFile
dark_contigs_file=args.unknownContigFile
runid=args.runid
bioproject=args.bioproject
output=args.outputFile

#Function to check whether input is fasta
def is_fasta(filename):
    with open(filename,'r') as handle:
        fasta = SeqIO.parse(handle,'fasta')
        return any(fasta)

def get_fasta_len_nuc(fastafile):
    total_bases=0
    for record in SeqIO.parse(fastafile, "fasta"):
        #print(len(record.seq))
        total_bases=total_bases+len(record.seq)
        #print(f'Number of bases in {fastafile}: {total_bases}')
    return total_bases


input_fasta_files_list=[rawcontigfile,long_contigs_file,grey_contigs_file,dark_contigs_file]

for fasta_file in input_fasta_files_list:
    if (is_fasta(fasta_file)):
        pass
    else:
        print('Input a valid FASTA file \n')

#for fasta_file in input_fasta_files_list:
#    get_fasta_len_nuc(fasta_file)

raw_contigs_bases=get_fasta_len_nuc(rawcontigfile)
long_contigs_bases=get_fasta_len_nuc(long_contigs_file)
grey_contigs_bases=get_fasta_len_nuc(grey_contigs_file)
contigs_unknown_origin_bases=get_fasta_len_nuc(dark_contigs_file)

if float(contigs_unknown_origin_bases) == 0 or float(long_contigs_bases) == 0:
    unknown_contig_bases_prop=0
else:
    unknown_contig_bases_prop=round((float(contigs_unknown_origin_bases)*100)/float(long_contigs_bases),2)

outfile=open(output,'w')
outfile.write('BioProject'+','
              +'SRARunID'+','
              +'TotalBases'+','
              +'LongBases(>=300_contigs)'+','
              +'PartiallyKnownBases(>=300_contigs)'+','
              +'BasesWithUnknownOrigin'+','
              +'UnknownBasesProportion'+'\n')

outfile.write(bioproject+','
          +runid+','
          +str(raw_contigs_bases)+','
          +str(long_contigs_bases)+','
          +str(grey_contigs_bases)+','
          +str(contigs_unknown_origin_bases)+','
          +str(unknown_contig_bases_prop)+'\n')
outfile.close()
