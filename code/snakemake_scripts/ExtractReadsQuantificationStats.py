#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon 16 Dec 12:24:36 GMT 2019

@author: sejmodha
"""
import argparse, subprocess, gzip

parser = argparse.ArgumentParser(description='This script parses raw FASTQ files and mapped .bam files to generate reads statistics')
parser.add_argument('-f1','--rawReadsFile', help='Input a FASTQ file',required=True)
parser.add_argument('-f2','--rawReadsFile2', help='Input a FASTQ file',required=False)
parser.add_argument('-u','--humanBAMFile',help='Input a valid BAM file of reads mapped to human genome', required=True)
parser.add_argument('-c','--contigsBAMFile',help='Input a valid BAM file of reads mapped to assembled contigs',required=True)
parser.add_argument('-r','--runid',help='ENA/SRA run id',required=True)
parser.add_argument('-b','--bioproject',help='BioProject number',required=True)
parser.add_argument('-o','--outputFile',help='Specify a name for the output file',required=True)
parser.add_argument('-t','--threads',help='Number of CPUs to use',required=False)
args=parser.parse_args()

rawreadsfile1=args.rawReadsFile
rawreadsfile2=args.rawReadsFile2
humanbamfile=args.humanBAMFile
contigsbamfile=args.contigsBAMFile
runid=args.runid
bioproject=args.bioproject
output=args.outputFile
threads=args.threads

if threads is None:
    threads=str(2)

#print(threads)

#function to count number of records in fastq file(s) using readlength.sh from bbtools
def count_FASTQ_records(fastqfile1,fastqfile2=None):
    if fastqfile2 is None:
        len_out=subprocess.check_output('readlength.sh -qin=64 -in='+fastqfile1,shell=True)
    else:
        len_out=subprocess.check_output('readlength.sh -in='+fastqfile1+' -in2='+fastqfile2,shell=True)
    rawreads=len_out.decode(encoding='utf-8').split('\n')[0].split('\t')[1]
    rawbases=len_out.decode(encoding='utf-8').split('\n')[1].split('\t')[1]
    avglen=len_out.decode(encoding='utf-8').split('\n')[4].split('\t')[1]
    return [rawreads,rawbases,avglen]


#function to extract reads data from BAM file
def get_BAM_reads_stats(bamfile):
    stats=subprocess.check_output('samtools stats -@ '+threads+' '+bamfile+'|grep ^SN | cut -f 2-',shell=True)
    stats=stats.decode(encoding='utf-8').split('\n')
    total_reads=stats[0].split('\t')[1]
    total_bases=stats[15].split('\t')[1]
    mapped_reads=stats[6].split('\t')[1]
    unmapped_reads=stats[8].split('\t')[1]
    mapped_bases=stats[18].split('\t')[1]
    unmapped_bases=int(total_bases)-int(mapped_bases)
    avg_len_reads=stats[24].split('\t')[1]
    return [total_reads,total_bases,mapped_reads,unmapped_reads,mapped_bases,unmapped_bases,avg_len_reads]

raw_reads_stats=count_FASTQ_records(rawreadsfile1,rawreadsfile2)
human_bam_stats=get_BAM_reads_stats(humanbamfile)
contig_bam_stats=get_BAM_reads_stats(contigsbamfile)
#
# print(raw_stats)
# print(human_bam_stats)
# print(contig_bam_stats)

outfile=open(output,'w')
outfile.write('BioProject'+','
              +'SRARunID'+','
              +'TotalRawReads_SRA'+','
              +'TotalRawBases_SRA'+','
              +'AverageReadLen_SRA'+','
              +'TotalReads_Human'+','
              +'TotalBases_Human'+','
              +'MappedReads_Human'+','
              +'UnmappedReads_Human'+','
              +'MappedBases_Human'+','
              +'UnmappedBases_Human'+','
              +'AverageReadLen_Human'+','
              +'TotalReads_Contigs'+','
              +'TotalBases_Contigs'+','
              +'MappedReads_Contigs'+','
              +'UnmappedReads_Contigs'+','
              +'MappedBases_Contigs'+','
              +'UnmappedBases_Contigs'+','
              +'AverageReadLen_Contigs'
              +'\n')

outfile.write(bioproject+','
          +runid+','
          +str(raw_reads_stats[0])+','
          +str(raw_reads_stats[1])+','
          +str(raw_reads_stats[2])+','
          +str(human_bam_stats[0])+','
          +str(human_bam_stats[1])+','
          +str(human_bam_stats[2])+','
          +str(human_bam_stats[3])+','
          +str(human_bam_stats[4])+','
          +str(human_bam_stats[5])+','
          +str(human_bam_stats[6])+','
          +str(contig_bam_stats[0])+','
          +str(contig_bam_stats[1])+','
          +str(contig_bam_stats[2])+','
          +str(contig_bam_stats[3])+','
          +str(contig_bam_stats[4])+','
          +str(contig_bam_stats[5])+','
          +str(contig_bam_stats[6])
          +'\n')

outfile.close()
