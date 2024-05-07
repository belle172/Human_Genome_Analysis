# -*- coding: utf-8 -*-
''' 
Created Dec 05 2023 
Python 3.10 
Jasper Bellefeuille - belle172@umn.edu 
Repository: Human_Genome_Analysis/vcf_genome_reader.py 

This script processes a vcf genome file into 

Documentation of the vcf genome file type: 
    support.illumina.com/help/BaseSpace_App_TumorNormal_help/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_gVCF.htm

Columns: CHROM - POS (position in the reference chromosome) 
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
# ref_allele = 0 
for line in vcf_file: 

    if line[0:2] == '##': 
        info_file.write(line) 
        # variants_file.write(line) 

    elif line[0:2] == '#C': 
        variants_file.write(line) 
        headers = line.strip().split('\t') 
        form_col = headers.index('FORMAT') 

    else: # not a header line 
        data = line.strip().split('\t') 

        # non-data entries to be filtered: GT = ./. means genotype could not be determined, and  
        if data[8][0:5] == 'GT:DP' and data[9][0:5] == './.:0':         # DP = 0 means no reads
            non_data += 1 

        ## Uncomment this block if you don't want to include reference alleles in the genotype 
        ## genotype is homozygous reference allele: GP = 0/0, with probability GP = 1 
        # elif data[8][0:5] == 'GT:GP' and data[9][0:5] == '0/0:1': 
        #     ref_allele += 1 

        else: # variant! 
            variants_file.write(line) 
            if data[8][0:5] == 'GT:DP': 
                reads = data[9].split(':') 
                chroms[data[0]] += int(reads[1]) 

vcf_file.close() 
info_file.close() 
variants_file.close() 

# Write a smaller file with a few of the entries in the genome file 
variants_file = open('genome_variants.txt') 
small_file = open('vcf_start.txt', 'w', encoding='utf-8') 

for line in variants_file: 
    data = line.strip().split('\t') 
    if data[7] != 'LOWQ' and data[9][0:3] != '0/0' and data[9][0:3] != './.': 
        small_file.write(line) 

small_file.close() 
variants_file.close() 

