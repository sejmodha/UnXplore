import os

configfile:
    #Config file for snakemake workflow
    "/path/to/sample_config.json"
    
# A Snakemake regular expression matching the forward mate FASTQ files.
DIRS, SAMPLES, = glob_wildcards(config['data']+"/{dir}/{sample}_1.fastq.gz")

#print(DIRS)
#Subset the DIRS list to remove any directories that may have been generated in the previous snakemake run
DIRS = [x for x in DIRS if "spades" not in x]
print(DIRS)
print(SAMPLES)
#Subset the SAMPLES list to remove any files that may have been generated in the previous snakemake run
SAMPLES = [x for x in SAMPLES if "GRCh38" not in x]
#print(SAMPLES)

rule all:
    input:
        #expand(config['data']+'/{dir}/{sample}_trim_1.fastq',zip,dir=DIRS,sample=SAMPLES),
        #expand(config['data']+'/{dir}/{sample}_GRCh38mapped.bam',zip,dir=DIRS,sample=SAMPLES)
        #expand(config['data']+'/{dir}/{sample}_GRCh38_unmapped_1.fq',zip,dir=DIRS,sample=SAMPLES)
        #expand(directory(config['data']+'/{dir}/{sample}_spades'),zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_contigs_mapped.bam',zip, dir=DIRS,sample=SAMPLES),
        # expand(config['data']+'/{dir}/{sample}_spades/contigs.fasta',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_diamond.m8',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_blastx_unaligned.fasta',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_blastn.m7',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_contigs_mapped.bam',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_dark_contigs.fasta',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_bbmap_pileup.txt',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_long_contigs_bbmap_pileup.txt',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_dark_contigs_bbmap_pileup.txt',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_contig_quantification.csv',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_diamond_partially_known.m8',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_bases_quantification.csv',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_partially_known_diam_pileup_stats.txt',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_partially_known_LCA.tsv',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_all_diam_pileup_stats.txt',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_all_diam_known_LCA.tsv',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_viral_lca_pileup.txt',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_viral_contigs.fa',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_viral_contigs_blastn.m6',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_viral_contigs_blastn_LCA.tsv',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_viral_contigs_blastn_LCA_metadata.csv',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_reads_stats.csv',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_diamond_known.m8',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_known_contigs.fasta',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_known_contig_names.txt',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_all_blastn_classified_LCA.tsv',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_diam_KNOWN_viral_contigs_metadata.csv',zip,dir=DIRS,sample=SAMPLES),
        expand(config['data']+'/{dir}/{sample}_diam_KNOWN_viral_contigs.fa',zip,dir=DIRS,sample=SAMPLES)

rule trimming:
    input:
        r1=config['data']+'/{dir}/{sample}_1.fastq.gz',
        r2=config['data']+'/{dir}/{sample}_2.fastq.gz'
    output:
        r1 = temp(config['data']+'/{dir}/{sample}_trim_1.fq.gz'),
        r2 = temp(config['data']+'/{dir}/{sample}_trim_2.fq.gz')
    params:
        adapters=config['adaptors'],
        phix_ill=config['phix_ill'],
        phix_adapt=config['phix_adapt'],
        lambda_seq=config['lambda'],
        pjet1=config['pjet1'],
        small=config['small']
    threads: 12
    message: ''' --- Trimming  --- '''
    conda:
        "phase1.yml"
    shell:
        """
        bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2}\
         ref={params.adapters},{params.phix_ill},{params.phix_adapt},{params.lambda_seq},{params.pjet1},{params.small}\
         ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=20 tpe tbo t={threads} -Xmx32g
        """

rule human_unmapped:
    input:
        r1 = config['data']+'/{dir}/{sample}_trim_1.fq.gz',
        r2 = config['data']+'/{dir}/{sample}_trim_2.fq.gz'
    output:
        bam = config['data']+'/{dir}/{sample}_GRCh38mapped.bam',
        r1 = temp(config['data']+'/{dir}/{sample}_GRCh38_unmapped_1.fastq'),
        r2 = temp(config['data']+'/{dir}/{sample}_GRCh38_unmapped_2.fastq')
    params:
        ref=config['genome']
    threads: 6
    conda:
        "phase1.yml"
    shell:
        """
        bwa mem {params.ref} {input} -t {threads}|samtools view -b -h -@ {threads}| samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        samtools fastq -1 {output.r1} -2 {output.r2} -n -f 4 -0 /dev/null -s /dev/null {output.bam} -@ {threads}
        """

rule bbmap_diginorm:
    input:
        r1 = config['data']+'/{dir}/{sample}_GRCh38_unmapped_1.fastq',
        r2 = config['data']+'/{dir}/{sample}_GRCh38_unmapped_2.fastq'
    output:
        r1 = temp(config['data']+'/{dir}/{sample}_GRCh38_unmapped_diginorm_1.fq.gz'),
        r2 = temp(config['data']+'/{dir}/{sample}_GRCh38_unmapped_diginorm_2.fq.gz')
    params:
        target=100,
        hist=config['data']+'/{dir}/{sample}_bbmap_in_hist.txt',
        histout=config['data']+'/{dir}/{sample}_bbmap_out_hist.txt',
        mindepth=3
    threads: 6
    conda:
        "phase1.yml"
    shell:
        """
        bbnorm.sh in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} target={params.target} hist={params.hist} histout={params.histout} t={threads} mindepth={params.mindepth} -Xmx32g
        """

rule spades_meta_assemble:
    input:
        r1 = config['data']+'/{dir}/{sample}_GRCh38_unmapped_diginorm_1.fq.gz',
        r2 = config['data']+'/{dir}/{sample}_GRCh38_unmapped_diginorm_2.fq.gz'
    output:
        contigs = config['data']+'/{dir}/{sample}_spades/contigs.fasta'
    params:
        outdir = directory(config['data']+'/{dir}/{sample}_spades'),
        ram = 128
    #threads: 8
    conda:
        "phase1.yml"
    shell:
        """
        spades.py --meta -1 {input.r1} -2 {input.r2} -o {params.outdir} -m {params.ram}
        """

rule bwa_contig_map:
    input:
        ref = config['data']+'/{dir}/{sample}_spades/contigs.fasta',
        r1 = config['data']+'/{dir}/{sample}_GRCh38_unmapped_1.fastq',
        r2 = config['data']+'/{dir}/{sample}_GRCh38_unmapped_2.fastq'
    output:
        config['data']+'/{dir}/{sample}_contigs_mapped.bam'
    params:
        #ref = config['data']+'/{dir}/{sample}_spades/contigs.fasta',
        tempref = config['data']+'/{dir}/{sample}_contigs.fa'
    conda:
        "phase1.yml"
    shell:
        """
        cp {input.ref} {params.tempref}
        bwa index {params.tempref}
        bwa mem {params.tempref} {input.r1} {input.r2} |samtools view -b -h | samtools sort -o {output} -
        samtools index {output}
        rm -f {params.tempref}*
        """

rule extract_long_contigs:
    input:
        contigs = config['data']+'/{dir}/{sample}_spades/contigs.fasta'
    output:
        config['data']+'/{dir}/{sample}_long_contigs.fasta'
    params:
        filter_fasta=config['customscripts']+'FilterFasta.py',
        length=300,
        prefix='{sample}'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.filter_fasta} -i {input} -o {output} -l {params.length} -p {params.prefix}
        """

rule diamond_blastx:
    input:
        contigs = config['data']+'/{dir}/{sample}_long_contigs.fasta'
    output:
        config['data']+'/{dir}/{sample}_diamond.m8'
    params:
        db = config['diamonddb']
    threads: 12
    conda:
        "phase1.yml"
    shell:
        """
        if [ -s {input} ];
        then
            diamond blastx -d {params.db} -p {threads} -q {input} -o {output} --unal 1 -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe staxids stitle
        else
            echo "{input} is empty - generating a pseudo file"
            touch {output} && chmod 755 {output}
        fi
        """

rule extract_diamond_unaligned:
    input:
        config['data']+'/{dir}/{sample}_diamond.m8'
    output:
        contigs = config['data']+'/{dir}/{sample}_blastx_unaligned.fasta'
    params:
        parse_blast_output = config['customscripts']+'ParseDIAMONDOutput.py',
        fasta = config['data']+'/{dir}/{sample}_long_contigs.fasta'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.parse_blast_output} -i {input} -f {params.fasta} -o {output}
        """

rule blastn:
    input:
        config['data']+'/{dir}/{sample}_blastx_unaligned.fasta'
    output:
        config['data']+'/{dir}/{sample}_blastn.m7'
    threads: 12
    params:
        evalue = 0.001,
        db = config['blastndb']
    conda:
        "phase1.yml"
    shell:
        """
        blastn -query {input} -db {params.db} -out {output} -evalue {params.evalue} -num_threads {threads} -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe staxids stitle"
        """

rule extract_blastn_unaligned:
    input:
        config['data']+'/{dir}/{sample}_blastn.m7'
    output:
        contigs = config['data']+'/{dir}/{sample}_dark_contigs.fasta'
    params:
        parse_blast_output = config['customscripts']+'ParseBLASTNOutput.py',
        fasta = config['data']+'/{dir}/{sample}_long_contigs.fasta'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.parse_blast_output} -i {input} -f {params.fasta} -o {output}
        """

rule get_assembly_stats_contigs:
    input:
        config['data']+'/{dir}/{sample}_contigs_mapped.bam'
    output:
        config['data']+'/{dir}/{sample}_bbmap_pileup.txt'
    params:
        contigs=config['data']+'/{dir}/{sample}_spades/contigs.fasta'
    conda:
        "phase1.yml"
    shell:
        """
        pileup.sh in={input} ref={params.contigs} out={output} 32bit=t overwrite=t -Xmx16g
        """

rule subset_pileup_stats:
    input:
        in_pileup=config['data']+'/{dir}/{sample}_bbmap_pileup.txt',
        long_contigs=config['data']+'/{dir}/{sample}_long_contigs.fasta',
        dark_contigs=config['data']+'/{dir}/{sample}_dark_contigs.fasta'
    params:
        get_pileup_subset=config['customscripts']+'GetPileupSubset.py',
        runid='{sample}',
        bioproject='{dir}'
    output:
        long_contigs_pileup=config['data']+'/{dir}/{sample}_long_contigs_bbmap_pileup.txt',
        dark_contigs_pileup=config['data']+'/{dir}/{sample}_dark_contigs_bbmap_pileup.txt'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.get_pileup_subset} -i {input.in_pileup} -f {input.long_contigs} -o {output.long_contigs_pileup} -r {params.runid} -b {params.bioproject}
        python3 {params.get_pileup_subset} -i {input.in_pileup} -f {input.dark_contigs} -o {output.dark_contigs_pileup} -r {params.runid} -b {params.bioproject}
        """

rule get_contigs_quantification:
    input:
        all_contigs=config['data']+'/{dir}/{sample}_spades/contigs.fasta',
        long_contigs=config['data']+'/{dir}/{sample}_long_contigs.fasta',
        grey_contigs=config['data']+'/{dir}/{sample}_grey_contigs.fasta',
        dark_contigs=config['data']+'/{dir}/{sample}_dark_contigs.fasta',
        diamond_out=config['data']+'/{dir}/{sample}_diamond.m8',
        blastn_out=config['data']+'/{dir}/{sample}_blastn.m7'
    params:
        get_contigs_quantification=config['customscripts']+'GetContigQuantificationStats.py',
        runid='{sample}',
        bioproject='{dir}'
    output:
        config['data']+'/{dir}/{sample}_contig_quantification.csv'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.get_contigs_quantification} -c {input.all_contigs} -l {input.long_contigs} -g {input.grey_contigs} -u {input.dark_contigs} -n {input.blastn_out} -d {input.diamond_out} -r {params.runid} -b {params.bioproject} -o {output}
        """

rule extract_partial_matches:
    input:
        diamond_file=config['data']+'/{dir}/{sample}_diamond.m8',
        fasta_file=config['data']+'/{dir}/{sample}_long_contigs.fasta'
    params:
        extract_partially_known=config['customscripts']+'ExtractPartiallyKnownSeq.py'
    output:
        diam_subset=config['data']+'/{dir}/{sample}_diamond_partially_known.m8',
        grey_fasta=config['data']+'/{dir}/{sample}_grey_contigs.fasta'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.extract_partially_known} -i {input.diamond_file} -f {input.fasta_file} -o {output.grey_fasta} -s {output.diam_subset}
        """

rule get_base_quantification:
    input:
        all_contigs=config['data']+'/{dir}/{sample}_spades/contigs.fasta',
        long_contigs=config['data']+'/{dir}/{sample}_long_contigs.fasta',
        dark_contigs=config['data']+'/{dir}/{sample}_dark_contigs.fasta',
        grey_contigs=config['data']+'/{dir}/{sample}_grey_contigs.fasta'
    params:
        get_bases_quantification=config['customscripts']+'GetBasesQuantification.py',
        runid='{sample}',
        bioproject='{dir}'
    output:
        config['data']+'/{dir}/{sample}_bases_quantification.csv'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.get_bases_quantification} -c {input.all_contigs} -l {input.long_contigs} -g {input.grey_contigs} -u {input.dark_contigs} -r {params.runid} -b {params.bioproject} -o {output}
        """

rule get_partially_known_stats:
    input:
        diam_subset=config['data']+'/{dir}/{sample}_diamond_partially_known.m8',
        long_contigs_pileup=config['data']+'/{dir}/{sample}_long_contigs_bbmap_pileup.txt'
    params:
        get_partially_known_stats=config['customscripts']+'ExtractSeqStats.py'
    output:
        partially_known_stats=config['data']+'/{dir}/{sample}_partially_known_diam_pileup_stats.txt'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.get_partially_known_stats} -d {input.diam_subset} -p {input.long_contigs_pileup} -o {output}
        """

rule partially_known_LCA:
    input:
        #partially_known_stats=config['data']+'/{dir}/{sample}_partially_known_diam_pileup_stats.txt'
        partially_known_stats=config['data']+'/{dir}/{sample}_diamond_partially_known.m8'
    params:
        get_partially_known_LCA=config['customscripts']+'ExtractLCA.py'
    output:
        partially_known_LCA_file=config['data']+'/{dir}/{sample}_partially_known_LCA.tsv'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.get_partially_known_LCA} -i {input.partially_known_stats} -o {output}
        """

rule joinDIAMONDandPileup:
    input:
        diam_file=config['data']+'/{dir}/{sample}_diamond.m8',
        pileup_file=config['data']+'/{dir}/{sample}_long_contigs_bbmap_pileup.txt'
    params:
        get_all_known_stats=config['customscripts']+'ExtractSeqStats.py'
    output:
        config['data']+'/{dir}/{sample}_all_diam_pileup_stats.txt'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.get_all_known_stats} -d {input.diam_file} -p {input.pileup_file} -o {output}
        """

rule all_LCA:
    input:
        #diam_known_stats=config['data']+'/{dir}/{sample}_all_diam_pileup_stats.txt'
        diam_known_stats=config['data']+'/{dir}/{sample}_diamond.m8'
    params:
        get_LCA=config['customscripts']+'ExtractLCA.py'
    output:
        config['data']+'/{dir}/{sample}_all_diam_known_LCA.tsv'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.get_LCA} -i {input.diam_known_stats} -o {output}
        """

rule extract_viral_hits:
    input:
        lca_file=config['data']+'/{dir}/{sample}_all_diam_known_LCA.tsv',
        pileup_file=config['data']+'/{dir}/{sample}_long_contigs_bbmap_pileup.txt',
        long_contigs=config['data']+'/{dir}/{sample}_long_contigs.fasta'
    params:
        extract_viral=config['customscripts']+'ExtractViralHits.py'
    output:
        outviral=config['data']+'/{dir}/{sample}_viral_lca_pileup.txt',
        viralcontigs=config['data']+'/{dir}/{sample}_viral_contigs.fa'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.extract_viral} -p {input.pileup_file} -f {input.long_contigs} -l {input.lca_file} -o {output.outviral} -s {output.viralcontigs}
        """

rule blastn_viral_contigs:
    input:
        config['data']+'/{dir}/{sample}_viral_contigs.fa'
    output:
        config['data']+'/{dir}/{sample}_viral_contigs_blastn.m6'
    threads: 12
    params:
        evalue = 0.001,
        db = config['blastndb']
    conda:
        "phase1.yml"
    shell:
        """
        blastn -query {input} -db {params.db} -out {output} -evalue {params.evalue} -num_threads {threads} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe staxids stitle qcovs"
        """

rule lca_blastn_viral_contigs:
    input:
        config['data']+'/{dir}/{sample}_viral_contigs_blastn.m6'
    output:
        config['data']+'/{dir}/{sample}_viral_contigs_blastn_LCA.tsv'
    params:
        getBLASTNViralLCA=config['customscripts']+'ExtractLCABLASTM6.py'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.getBLASTNViralLCA} -i {input} -o {output}
        """

rule lca_qcov_join:
    input:
        virallcafile=config['data']+'/{dir}/{sample}_viral_contigs_blastn_LCA.tsv',
        viralblastnfile=config['data']+'/{dir}/{sample}_viral_contigs_blastn.m6'
    output:
        config['data']+'/{dir}/{sample}_viral_contigs_blastn_LCA_metadata.csv'
    params:
        joinLCABLASTN=config['customscripts']+'JoinLCAAndQCovBLASTN.py',
        pileupfile=config['data']+'/{dir}/{sample}_long_contigs_bbmap_pileup.txt',
        fastafile=config['data']+'/{dir}/{sample}_viral_contigs.fa'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.joinLCABLASTN} -l {input.virallcafile} -b {input.viralblastnfile} -o {output} -p {params.pileupfile}  -f {params.fastafile}
        """

rule extract_reads_stats:
    input:
        fastq1=config['data']+'/{dir}/{sample}_1.fastq.gz',
        fastq2=config['data']+'/{dir}/{sample}_2.fastq.gz',
        humanbam=config['data']+'/{dir}/{sample}_GRCh38mapped.bam',
        contigbam=config['data']+'/{dir}/{sample}_contigs_mapped.bam'
    output:
        config['data']+'/{dir}/{sample}_reads_stats.csv'
    threads: 4
    params:
        get_reads_stats=config['customscripts']+'ExtractReadsQuantificationStats.py',
        runid='{sample}',
        bioproject='{dir}'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.get_reads_stats} -f1 {input.fastq1} -f2 {input.fastq2} -u {input.humanbam} -c {input.contigbam} -r {params.runid} -b {params.bioproject} -o {output} -t {threads}
        """

rule extract_known_matches:
    input:
        diamond_file=config['data']+'/{dir}/{sample}_diamond.m8',
        fasta_file=config['data']+'/{dir}/{sample}_long_contigs.fasta'
    params:
        extract_known=config['customscripts']+'ExtractKnownSeq.py'
    output:
        diam_known=config['data']+'/{dir}/{sample}_diamond_known.m8',
        known_fasta=config['data']+'/{dir}/{sample}_known_contigs.fasta',
        known_contig_names=config['data']+'/{dir}/{sample}_known_contig_names.txt'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.extract_known} -i {input.diamond_file} -f {input.fasta_file} -o {output.known_fasta} -s {output.diam_known} -n {output.known_contig_names}
        """

rule get_LCA_for_blastn:
    input:
        blastn_file=config['data']+'/{dir}/{sample}_blastn.m7'
    params:
        get_LCA=config['customscripts']+'ExtractLCA.py'
    output:
        config['data']+'/{dir}/{sample}_all_blastn_classified_LCA.tsv'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.get_LCA} -i <(grep -v '^#' {input.blastn_file}) -o {output}
        """

rule extract_known_viral_diamond:
    input:
        known_contig_names=config['data']+'/{dir}/{sample}_known_contig_names.txt',
        diam_classified_stats=config['data']+'/{dir}/{sample}_all_diam_pileup_stats.txt',
        diam_classified_lca=config['data']+'/{dir}/{sample}_all_diam_known_LCA.tsv',
        long_contigs_fasta=config['data']+'/{dir}/{sample}_long_contigs.fasta'
    params:
        get_viral_hits=config['customscripts']+'ExtractViralHitsDIAMEnhanced.py'
    output:
        known_viral_contigs_diam_lca=config['data']+'/{dir}/{sample}_diam_KNOWN_viral_contigs_metadata.csv',
        known_viral_contigs_fasta=config['data']+'/{dir}/{sample}_diam_KNOWN_viral_contigs.fa'
    conda:
        "phase1.yml"
    shell:
        """
        python3 {params.get_viral_hits} -n {input.known_contig_names} -f {input.long_contigs_fasta} -l {input.diam_classified_lca} -d {input.diam_classified_stats} -o {output.known_viral_contigs_diam_lca} -s {output.known_viral_contigs_fasta}
        """
onsuccess:
    print("Workflow finished, no error")

