#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon 23 Sep 11:24:38 BST 2019

@author: sejmodha
"""


import pandas as pd
import argparse
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()


parser = argparse.ArgumentParser(description='This script parses a tabualr BLAST output file and generates an output with LCA details')
parser.add_argument('-i','--infile', help='Input a valid file with following headers (qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qframe,staxids,stitle) ',required=True)
parser.add_argument('-o','--outfile',help='Output file name', required=True)

args=parser.parse_args()

inputfile = args.infile
outputfile = args.outfile


def generate_seqid_taxid_dict(inputfile):
    """A function that groups the input csv file according to subject taxonomy ID"""
    diam_header='qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qframe,staxids,stitle'.split(',')
    df = pd.read_csv(inputfile,sep='\t',low_memory=False,header=None,names=diam_header)
    #print(df.head())
    #Remove rows where taxonomy ID is unavailable
    notaxid_df = df[df['staxids'].isna()]
    df = df[df['staxids'].notna()]
    grouped_dict = df.groupby('qseqid')['staxids'].apply(list).to_dict()
    #print(grouped_dict)
    return grouped_dict

def get_taxonomy_id_list(inputdict):
    """A function that takes a dict argument and reformats the values in the list and saves it to a new dict"""
    formatted_taxa_dict = {}
    for seqid,taxids_list in inputdict.items():
        #if seqid == 'ERR537005_NODE_10046_length_3186_cov_4.572022':
        new_taxa_list = []
        for taxid in taxids_list:
            if ';' in str(taxid):
                #print(taxid)
                #new_taxa_list.extend(taxid.split(';'))
                new_taxa_list.append(taxid.split(';')[0])
            else:
                new_taxa_list.append(taxid)
        formatted_taxa_dict[seqid]=list(set(new_taxa_list))
    #print(formatted_taxa_dict)
    return formatted_taxa_dict

def get_ancestors(taxid):
    """A function to get lineage of a given taxonomy ID"""
    lineage = []
    try:
        lineage = ncbi.get_lineage(taxid)
    except ValueError:
        lineage = None
        pass
    return lineage

def get_lca(taxon1,taxon2):
    """Function to get LCA between two taxonomy IDs"""
    if taxon1 is not None and taxon2 is not None:
        viralhits = 0
        #print(f'getting LCA of {taxon1} and {taxon2}')
        ancestors_1 = get_ancestors(taxon1)[::-1]
        ancestors_2 = get_ancestors(taxon2)[::-1]
        for ancestor in ancestors_1:
            #print(ancestor)
            if ancestor in ancestors_2:
                return ancestor
    else:
        print(taxon1, taxon2)

def get_lca_list(taxalist):
    """A function to get LCA of a given list of taxonomy IDs"""
    #print(f'Taxalist is : {taxalist}')
    taxon1 = int(taxalist.pop())
    while len(taxalist) > 0:
        taxon2 = int(taxalist.pop())
        lca = get_lca(taxon1, taxon2)
        #print(f'LCA of {taxon1} and {taxon2} is {lca}')
        taxon1 = lca
    return taxon1

def get_viral_prop(infile,outfile):
    """Function to process inputfile with partial matches and produces
    an output file with proportion of viral hits for each contig ID"""
    output=open(outfile,'w')
    output.write('ContigID\tLCA_TaxonID\tLCA_TaxonName\tNumberOfUniqueTaxIDs\tNumberOfViralHits\tProportionOfViralHits\n')
    for seqid,taxaid_list in get_taxonomy_id_list(generate_seqid_taxid_dict(infile)).items():
        #taxa_lineage_list = list(map(lambda x: ncbi.get_lineage(int(x)), taxaid_list))
        taxa_lineage_list = list(map(lambda x: get_ancestors(int(x)), taxaid_list))
        print(taxa_lineage_list)
        #print(seqid,taxa_lineage_list,len(taxa_lineage_list))
        lca = {}
        virus_count = 0
        viral_prop = 0
        if all(x is None for x in taxa_lineage_list):
            print(f'{seqid},{taxa_lineage_list} contains only NULL values, LCA cannot be computed')
        elif any(10239 == x[1] for x in taxa_lineage_list):
            virus_count = sum(x.count(10239) for x in taxa_lineage_list)
            viral_prop = round(virus_count*100/len(taxa_lineage_list),2)
            lca = ncbi.get_taxid_translator([get_lca_list(taxaid_list)])
            output.write(seqid+'\t'+str(list(lca.keys())[0])+'\t'+str(list(lca.values())[0])+'\t'+str(len(taxa_lineage_list))+'\t'+str(virus_count)+'\t'+str(viral_prop)+'\n')

        else:
            lca = ncbi.get_taxid_translator([get_lca_list(taxaid_list)])
            output.write(seqid+'\t'+str(list(lca.keys())[0])+'\t'+str(list(lca.values())[0])+'\t'+str(len(taxa_lineage_list))+'\t'+str(virus_count)+'\t'+str(viral_prop)+'\n')
    output.close()


get_viral_prop(inputfile,outputfile)
