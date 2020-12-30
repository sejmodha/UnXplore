#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon 14 Oct 15:36:23 BST 2019

@author: sejmodha
"""

import pandas as pd
import argparse
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

parser = argparse.ArgumentParser(description='This script parses a BLAST tabular output file generates an output with LCA details')
parser.add_argument('-i','--infile', help='Input a valid file in BLASTN format M6',required=True)
parser.add_argument('-o','--outfile',help='Output file name', required=True)
parser.add_argument('-c','--columns',help='Input column names',required=False,default='qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qframe,staxids,stitle,qcovs')


args=parser.parse_args()

inputfile = args.infile
outputfile = args.outfile
header=args.columns
print('Assuming the input file column order to be: '+header)

def generate_seqid_taxid_dict(inputfile):
    """A function that groups the input csv file according to subject taxonomy ID"""
    df=pd.read_csv(inputfile,sep='\t',low_memory=False,names=header.split(','))
    #print(df.head())
    #Remove rows where taxonomy ID is unavailable
    notaxid_df=df[df['staxids'].isna()]
    df=df[df['staxids'].notna()]
    grouped_dict = df.groupby('qseqid')['staxids'].apply(list).to_dict()
    #print(f'Length of grouped is {len(grouped)}')
    return grouped_dict

def get_taxonomy_id_list(inputdict):
    """A function that takes a dict argument and reformats the values in the list and saves it to a new dict"""
    formatted_taxa_dict={}
    for seqid,taxids_list in inputdict.items():
        #if seqid == 'ERR537005_NODE_10046_length_3186_cov_4.572022':
        new_taxa_list=[]
        for taxid in taxids_list:
            if ';' in str(taxid):
                #print(taxid)
                #new_taxa_list.extend(taxid.split(';'))
                new_taxa_list.append(taxid.split(';')[0])
            else:
                new_taxa_list.append(taxid)
        formatted_taxa_dict[seqid]=list(set(new_taxa_list))
    return formatted_taxa_dict

def get_ancestors(taxid):
    """A function to get lineage of a given taxonomy ID"""
    lineage=[]
    try:
        lineage=ncbi.get_lineage(taxid)
    except ValueError:
        pass
    return lineage

def get_lca(taxon1,taxon2):
    """Function to get LCA between two taxonomy IDs"""
    if taxon1 is not None and taxon2 is not None:
        viralhits=0
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

def GetFamily(taxon):
    lineage=[]
    try:
        lineage=ncbi.get_lineage(taxon)
    except ValueError:
        pass
    for tax in lineage:
        #print(tax)
        rank=ncbi.get_rank([tax])
        #print(rank)
        if 'family' in rank.values():
            family=ncbi.get_taxid_translator(rank.keys())
            #print(type(superkingdom.values()))
            return list(family.values())[0]

def GetSuperkingdom(taxon):
    lineage=[]
    try:
        lineage=ncbi.get_lineage(taxon)
    except ValueError:
        pass
    for tax in lineage:
        #print(tax)
        rank=ncbi.get_rank([tax])
        #print(rank)
        if 'superkingdom' in rank.values():
            superkingdom=ncbi.get_taxid_translator(rank.keys())
            #print(type(superkingdom.values()))
            return list(superkingdom.values())[0]


def get_viral_prop(infile,outfile):
    """Function to process inputfile with partial matches and produces
    an output file with proportion of viral hits for each contig ID"""
    output=open(outfile,'w')
    output.write('ContigID\tLCA_TaxonID\tLCA_TaxonName\tNumberOfUniqueTaxIDs\tNumberOfViralHits\tProportionOfViralHits\tLCA_Family\tLCA_Superkingdom\n')
    for seqid,taxaid_list in get_taxonomy_id_list(generate_seqid_taxid_dict(infile)).items():
        taxa_lineage_list=list(map(lambda x: ncbi.get_lineage(int(x)), taxaid_list))
        #print(seqid,taxa_lineage_list)
        lca={}
        virus_count=0
        viral_prop=0
        if any(10239 == x[1] for x in taxa_lineage_list):
            virus_count=sum(x.count(10239) for x in taxa_lineage_list)
            viral_prop=round(virus_count*100/len(taxa_lineage_list),2)
            lca=ncbi.get_taxid_translator([get_lca_list(taxaid_list)])
            family=GetFamily(list(lca.keys())[0])
            superkingdom=GetSuperkingdom(list(lca.keys())[0])
            #print(family, superkingdom)
            output.write(seqid+'\t'+str(list(lca.keys())[0])+'\t'+str(list(lca.values())[0])+'\t'+str(len(taxa_lineage_list))+'\t'+str(virus_count)+'\t'+str(viral_prop)+'\t'+str(family)+'\t'+str(superkingdom)+'\n')

        else:
            lca=ncbi.get_taxid_translator([get_lca_list(taxaid_list)])
            family=GetFamily(list(lca.keys())[0])
            superkingdom=GetSuperkingdom(list(lca.keys())[0])
            #print(family, superkingdom)
            output.write(seqid+'\t'+str(list(lca.keys())[0])+'\t'+str(list(lca.values())[0])+'\t'+str(len(taxa_lineage_list))+'\t'+str(virus_count)+'\t'+str(viral_prop)+'\t'+str(family)+'\t'+str(superkingdom)+'\n')
    output.close()


get_viral_prop(inputfile,outputfile)
