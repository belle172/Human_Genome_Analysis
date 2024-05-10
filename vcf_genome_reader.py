# -*- coding: utf-8 -*-
''' 
Python 3.10 
Created Dec 05 2023 
Jasper Bellefeuille - belle172@umn.edu 
Repository: Human_Genome_Analysis/vcf_genome_reader.py 

This script processes a vcf genome file to create a SNPome file usable for the Human Genome 
    Analysis done by phenotype_predictor.py. 

Assumptions for input data: 
    Files are not compressed 
    The vcf genome file uses hg38 as its reference genome 
    The vcf file does not provide the rsids of each polymorphism 
    The hg38 chromosomal position is provided for each line in the vcf file 

Documentation of the vcf genome file type: 
    support.illumina.com/help/BaseSpace_App_TumorNormal_help/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_gVCF.htm
    https://samtools.github.io/hts-specs/VCFv4.3.pdf

Columns: POS - position in the reference chromosome 
         REF ALT QUAL FILTER INFO FORMAT US-PU7C-[genome] 
''' 

# Assumes reference genome used is human genome 38, stated as assembly=hg38 in the header data 
chroms = {'chr1':0, 'chr2':0, 'chr3':0, 'chr4':0, 'chr5':0, 'chr6':0, 'chr7':0, 'chr8':0, 
          'chr9':0, 'chr10':0, 'chr11':0, 'chr12':0, 'chr13':0, 'chr14':0, 'chr15':0, 'chr16':0, 
          'chr17':0, 'chr18':0, 'chr19':0, 'chr20':0, 'chr21':0, 'chr22':0, 'chrX':0, 'chrY':0, 
          'chrM':0 } 

vcf_file = open('US-PU7C.vcf') 
info_file = open('vcf_metadata.txt', 'w', encoding='utf-8') 
variants_file = open('genome_variants.txt', 'w', encoding='utf-8') 

non_data = 0 
for line in vcf_file: 

    if line[0:2] == '##': # write the meta information header lines to a file 
        info_file.write(line) 
        # variants_file.write(line) # uncomment to include the headers 

    elif line[0:2] == '#C': # just include the column header line in the new genome file 
        headers = line.strip().split('\t') 
        form_col = headers.index('FORMAT') 

        skip = headers.index('ID') 
        s = '' 
        for header in headers: 
            if headers.index(header) not in skip: 
                s += header 
        variants_file.write(line) 

    # For each polymorphism line, only include if it has any reads 
    else: # not a header line 
        data = line.strip().split('\t') 

        # non-data entries to be filtered: DP = 0 means no reads 
        if data[8][0:5] == 'GT:DP' and data[9][0:5] == './.:0': 
            non_data += 1 

        # GT = . means a call cannot be made for that allele 
        # GQ != 0 means there is at least one alt allele call 
        # GQ = difference between lowest likelyhood genotype and second lowest likelyhood genotype 
        # −10log10 phred quality probability call is wrong conditioned on being variant 

        else: # variant! 

            # TODO Data cleanup: ID column has no information 
            # <*> is preferred over <NON_REF> 

            variants_file.write(line) 
            if data[8][0:5] == 'GT:DP': 
                reads = data[9].split(':') 
                chroms[data[0]] += int(reads[1]) 

vcf_file.close() 
info_file.close() 
variants_file.close() 

# =============================================================================
# Write a smaller file with a few of the entries in the genome file 
# =============================================================================
variants_file = open('genome_variants.txt') 
small_file = open('vcf_start.txt', 'w', encoding='utf-8') 

i=0 
for line in variants_file: 
    if i < 100000: 
        i += 1 # limit file size 
        data = line.strip().split('\t') 
        if data[6] != 'LOWQ' and data[9][0:3] != '0/0' and data[9][0:3] != './.': 
            small_file.write(line) 

small_file.close() 
variants_file.close() 

