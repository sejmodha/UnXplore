    ########################################################
    ########################################################
    ##                                                    ##
    ##       _   _      __  __      _                     ##
    ##      | | | |_ __ \ \/ /_ __ | | ___  _ __ ___      ##
    ##      | | | | '_ \ \  /| '_ \| |/ _ \| '__/ _ \     ##
    ##      | |_| | | | |/  \| |_) | | (_) | | |  __/     ##
    ##       \___/|_| |_/_/\_\ .__/|_|\___/|_|  \___|     ##
    ##                       |_|                          ##
    ##                                                    ##
    ########################################################


# UnXplore

[![DOI](https://zenodo.org/badge/324996408.svg)](https://zenodo.org/badge/latestdoi/324996408)


#### Exploring the unknown sequence matter embedded in human microbiome data

Mining human microbiomes for taxonomically unknown sequences.

<!--- See documentation available on [wiki](https://github.com/sejmodha/Unknown-Sequences/wiki) pages.  --->

Preprint: [Quantifying and cataloguing unknown sequences within human microbiomes](https://doi.org/10.1101/2021.01.22.427751) 

# Table of Contents
1. [Unknown sequence mining framework](#unknown-sequence-mining-framework)
2. [Framework application](#framework-application)


## Unknown sequence mining framework

In this study, we have developed a systematic framework to identify the so-called ‘dark’ sequence matter hidden within the publicly available human metagenomic datasets. This framework automatically cleans and filters short read sequence datasets before assembling them into larger contigs for taxonomic classification and investigation. This comprehensive analyses framework can identify taxonomically unknown contigs and can be applied to any human microbiome study.

A brief overview of the analyses steps is shown in the figure here:

![Framework](https://raw.githubusercontent.com/sejmodha/UnXplore/main/images/rulegraph.svg?token=ABO74WHS5J6JMZNRLIK254LABFPXI)

### Installation

#### Prerequisits:

You would require `miniconda` or `anaconda` installed on your machine.

1. Clone this repository on your local machine.
2. Create a new conda environment for program installation using `conda env create -f phase1.yml` - This should ensure that all programs required for the analyses are installed properly in this new environment.
3. Create a new environment for snakemake install using `conda env create -f snakemake.yml`
4. Download the following databases:
    - [BLAST `nt` databases](https://ftp.ncbi.nlm.nih.gov/blast/db/)
    - [BLAST `nr` databases in FASTA format](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/) - These FASTA sequences will need to be formatted into `diamond` compatible database 
    - `BWA` index of the human genome sequence 
5. Generate a `config` file in `JSON` format - A sample config file is available in this repository that can be modified with updated `PATH` for all variables included in the config file.
6. Select the single-end or paired-end analysis snakemake workflow to run the analysis pipeline.

**Note:** It is possible to install individual software/programs in their independent conda environments. If you have got that setup then you must:
- Export relevant conda environment in `.yaml` file and make that file available in the same folder as `.smk` file
and
- You must also change the following parameter in the `.smk` file for any rules that may use the software to
  ```
      conda:
        "myenv.yml"
  ```

A sample command to run the pipeline:

Dry run: `snakemake -j24 --use-conda -s PairedendAnalysis.smk -np`

Analysis: `snakemake -j24 --use-conda -s PairedendAnalysis.smk`


## Framework application

This framework has been applied to a large data set consisting of 963 metagenomic samples from 40 distinct studies that cover 10 different human microbiomes (ranging from gut and fecal, to oral and skin). This has enabled the identification of over half a million sequences that do not bear any sequence similarity to the sequences in the GenBank databases, so-called ‘dark’ sequences. We show that up to 25% of the assembled contigs within a microbiome are ‘dark’, particularly for the relatively unexplored human skin and oral microbiomes. 

A third of all unknown sequences (n=215,985) were shown to contain large predicted open reading frames (at least 100 amino acid long) and a small proportion of these open reading frames contained domain signatures confirming the presence of currently unidentified organisms. Moreover, a comprehensive clustering analysis has led to the identification of unknown sequences that were present across different human biomes (as well as from different samples/studies investigating the same human microbiome) indicating that we have discovered potentially widespread and as yet unclassified novel biological organisms within the human microbiome. Our analyses over the past 18 months, highlight that these dark sequences are being classified at a relatively low rate of 1.64% per month and our framework provides an approach to accelerate the rate by which these dark sequences can be identified and in turn be classified. This reiterates our initial hypothesis that these assembled unknown sequences are indeed of biological origin.

Analysis files generated from above are available in the results folder. 

### Results and output files:
**Note**: The first line of each file contains column names that are self-explanatory.  
- `UnknownContigsBasesResults.csv` - quantification of unknown contigs in each sample
- `ClusterAnalysis.csv` - `MMSeqs2` clustering results for unknown contigs
- `UnknownContigsAssemblyMetadata.csv` - assembly statistics for each unknown contig
- `FunctionalAnalysis.csv` - InterProScan analysis for predicted open reading frames generated from unknown contigs

### More details:
For further details about this study, see: https://doi.org/10.1101/2021.01.22.427751



