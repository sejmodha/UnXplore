#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 14:30:19 2019

@author: sejmodha
"""

import pandas as pd
import os, subprocess, argparse

parser = argparse.ArgumentParser(description='This script downloads reads metadata')
parser.add_argument('-i','--sraruninfo', help='Input SRA RunInfo file without header',required=True)
parser.add_argument('-l','--location',help='Specify output location', required=False,default=os.getcwd())

args = parser.parse_args()
#File download location
pwd=args.location

print('\nFiles will be saved in: '+pwd+'\n')

columns=['Run','ReleaseDate','LoadDate','spots','bases',
         'spots_with_mates','avgLength','size_MB','AssemblyName','download_path',
         'Experiment','LibraryName','LibraryStrategy','LibrarySelection','LibrarySource',
         'LibraryLayout','InsertSize','InsertDev','Platform','Model',
         'SRAStudy','BioProject','Study_Pubmed_id','ProjectID','Sample',
         'BioSample','SampleType','TaxID','ScientificName','SampleName',
         'g1k_pop_code','source','g1k_analysis_group','Subject_ID','Sex',
         'Disease','Tumor','Affection_Status','Analyte_Type','Histological_Type',
         'Body_Site','CenterName','Submission','dbgap_study_accession','Consent',
         'RunHash','ReadHash']

df = pd.read_csv(args.sraruninfo,header=None,names=columns)
#print(df.head(3))
#print(df['BioProject'])
#Save unique BioProject number in a list
bioproject=df.BioProject.unique().tolist()


#Iter through the list of the BioProject numbers
for i in range(len(bioproject)):
    os.chdir(pwd)
    print(bioproject[i])
    #Check if folders "PAIRED" and "SINGLE exists - if not then create them"
    if not os.path.exists('PAIRED'):
        os.makedirs('PAIRED')
    if not os.path.exists('SINGLE'):
        os.makedirs('SINGLE')
    #Save all RunID for a given BioProject in a list
    sra_run_id=list(df[df['BioProject'] == bioproject[i]]['Run'])
    #Identify the library layout for the sequencing
    library_layout=list(df[df['BioProject'] == bioproject[i]]['LibraryLayout'])[0]
    print(sra_run_id)
    print(library_layout)
    #Process this for the paired-end library
    for run_id in range(len(sra_run_id)):
        #Download paired-end data
        print("LOCATION IS: "+ os.getcwd())
        if 'PAIRED' in library_layout:
            os.chdir(pwd+'/PAIRED')
            #Check if the a folder with given BioProject number exists else create one
            if not os.path.exists(bioproject[i]):
                os.makedirs(bioproject[i])
                os.chdir(bioproject[i])
                print(f'prefetch + {sra_run_id[run_id]}')
                subprocess.call('prefetch  -X 999999999999999 '+str(sra_run_id[run_id]),shell=True)
                #subprocess.call('fastq-dump -A '+str(sra_run_id[run_id])+' --split-files',shell=True)
                subprocess.call('parallel-fastq-dump --sra-id '+str(sra_run_id[run_id])+' --threads 8 --split-files --gzip',shell=True)
            else:
                os.chdir(bioproject[i])
                subprocess.call('prefetch -X 999999999999999 '+str(sra_run_id[run_id]),shell=True)
                #subprocess.call('fastq-dump -A '+str(sra_run_id[run_id])+' --split-files',shell=True)
                subprocess.call('parallel-fastq-dump --sra-id '+str(sra_run_id[run_id])+' --threads 8 --split-files --gzip',shell=True)
        #Download single-end data
        else:
            os.chdir(pwd+'/SINGLE')
            #Check if the a folder with given BioProject number exists else create one
            if not os.path.exists(bioproject[i]):
                os.makedirs(bioproject[i])
                os.chdir(bioproject[i])
                subprocess.call('prefetch  -X 999999999999999 '+str(sra_run_id[run_id]),shell=True)
                #subprocess.call('fastq-dump -A '+str(sra_run_id[run_id]),shell=True)
                subprocess.call('parallel-fastq-dump --sra-id '+str(sra_run_id[run_id])+' --threads 8 --gzip',shell=True)
            else:
                os.chdir(bioproject[i])
                subprocess.call('prefetch  -X 999999999999999 '+str(sra_run_id[run_id]),shell=True)
                #subprocess.call('fastq-dump -A '+str(sra_run_id[run_id]),shell=True)
                subprocess.call('parallel-fastq-dump --sra-id '+str(sra_run_id[run_id])+' --threads 8 --gzip',shell=True)
