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
    Input vcf genome file uses human genome 38 (stated as assembly=hg38 in the header data) as 
    its reference genome and provides the hg38 chromosomal position for each line in the vcf file 

Documentation of the vcf genome file type: 
    support.illumina.com/help/BaseSpace_App_TumorNormal_help/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_gVCF.htm
    https://samtools.github.io/hts-specs/VCFv4.3.pdf

Columns: POS - position in the reference chromosome 
         REF ALT QUAL FILTER INFO FORMAT [genome] 
''' 

from file_processing import get_tsv_matrix 

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

    elif line[0:2] == '#C': # include the column header line in the new genome file 
        headers = line.strip().split('\t') 
        form_col = headers.index('FORMAT') 
        variants_file.write(line) 

        # TODO: use slim_matrix() to skip ID column 
        skip = headers.index('ID') 
        s = '' 
        for header in headers: 
            if headers.index(header) not in skip: 
                s += header 

    # For each polymorphism line, only include if it has any reads 
    else: # not a header line 
        data = line.strip().split('\t') 

        # non-data entries to be filtered: DP = 0 means no reads 
        # GT = . means a call cannot be made for that allele 
        if data[8][0:5] == 'GT:DP' and data[9][0:5] == './.:0': 
            non_data += 1 

        else: # variant! 
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
    data = line.strip().split('\t') 
    if data[0] != 'chr1': 
        if i < 100000: 
            i += 1 # limit file size 
            if data[6] != 'LOWQ' and data[9][0:3] != '0/0' and data[9][0:3] != './.': 
                small_file.write(line) 
        else: break 

small_file.close() 
variants_file.close() 


# =============================================================================
# Write SNPome file from a VCFv4.3 file of genome variants 
# =============================================================================
gwas_locs = get_tsv_matrix('gwas_snp_locations.txt') # chrom pos loc rsid 
chroms = {'1': 1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, 
          '12':12, '13':13, '14':14, '15':15, '16':16, '17':17, '18':18, '19':19, '20':20, 
          '21':21, '22':22, 'X':23, 'Y':24, 'M':25} 

variants_file = open('genome_variants.txt') 
variants_file.readline() 
snp_file = open('vcf_start.txt', 'w', encoding='utf-8') 
snp_file.write('#rsid\tchromosome\tposition\tgenotype\n') 
i_gwas = 0 

# Check if each variant location in the genome has a corresponding rsid 
for line in variants_file: 
    data = line.strip().split('\t') 
    chrom = chroms[data[0].lstrip('chr')] # chromosome of current line 
    try: 

        # catalog chrom < genome variant chrom 
        while int(gwas_locs[i_gwas][0]) < chrom: i_gwas += 1 

        # variant location is greater than catalog 
        while int(gwas_locs[i_gwas][1]) < int(data[1]): i_gwas += 1 

        if int(gwas_locs[i_gwas][1]) == int(data[1]): # write rsid chromosome position genotype 
            s = gwas_locs[i_gwas][3] + '\t' + data[0].lstrip('chr') + '\t' + data[1] + '\t' 
            gt = '' 

            if data[8][0:2] == 'GT' and data[9][0:3] != './.': # genotype is given 
                gt = data[9][0] + data[9][2] # change genotype format from './.' to '..' 

                # data validation 
                if data[9][0:3] not in ['0/0', '0/1', '0|1', '1|0', '1/1', '1|1', '1/2']: 
                    print(line) 

            # genotype isn't called (below 95% confidence) but genotype probability is given 
            elif data[8][0:5] == 'GT:GP': 
                gp = data[9].split(':')[1].split(',') 
                likely_gt = gp.index(max(gp)) 
                if likely_gt == 0: gt = '00' 
                if likely_gt == 1: gt = '10' 
                if likely_gt == 2: gt = '11' 

            # GQ = difference between lowest likelyhood genotype and second lowest likelyhood genotype 
            # −10log10 phred quality probability call is wrong conditioned on being variant 
            elif data[8][0:8] == 'GT:DP:GQ': # gene quality is given 
                gq = data[9].split(':')[2] 
                if gq == '0': gt = '00' # GQ == 0: all calls are reference. Set gt to ref/ref 

                # GQ != 0: means there is at least one alt allele call 
                else: 
                    s += line # there is at least one alt allele read 
                    print(line) 

            # TODO: parse all other cases 
            else: 
                s += line # genotype and probability aren't given 
                print(line) 

            for allele in gt: 
                if allele == '0': s += data[3] # reference allele, given in column 3 
                if allele == '1': s += data[4].split(',')[0] # first alt allele, given in column 4 
                if allele == '2': s += data[4].split(',')[1] # second alt allele 
            s += '\n' 

            snp_file.write(s) 
    except IndexError: break # if we've gone past the end of the catalog, stop comparing 

variants_file.close() 
snp_file.close() 

