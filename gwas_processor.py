# -*- coding: utf-8 -*- 
'''
Python 3.10 
Created Mar 30 2023 
Jasper Bellefeuille - belle172@umn.edu 
Repository: Human_Genome_Analysis/gwas_processor.py 

This script takes in the original GWAS catalog and outputs a smaller file containing a subset of 
    the columns and only SNPs with a specified allele, organized by rsid. 

GWAS catalog: downloaded from https://www.ebi.ac.uk/gwas/docs/file-downloads 
    All associations v1.0 
    Reference genome: assembly GRCh38.p14 and dbSNP Build 156 
        from FAQ page www.ebi.ac.uk/gwas/docs/faq#faq-E2 

Assumptions for input data: 
    The GWAS catalog is named gwas_catalog_v1-associations_e109.tsv 
    Catalog columns are the same as 03/30/2023 download 
''' 

import math # natural log function 
from file_processing import get_tsv_matrix, slim_matrix, rsid_matrix 

# GWAS catalog - each row is a genetic variant. file organized by ??, then each variant from the 
# same study is together 
gwas_matrix = get_tsv_matrix('gwas_catalog_v1-associations_e109.tsv') # read gwas catalog 

# Original gwas columns, notes, and indeces 
#   0:'DATE ADDED TO CATALOG', 1:'PUBMEDID' (also in LINK), 2:'AUTHOR', 3:'DATE', 4:'JOURNAL', 
#   5:'LINK' (link to pubmed entry), 6:'STUDY', 7:'TRAIT', 8:'SAMPLE SIZE', 
#   9:'REPLICATION SAMPLE SIZE', 10:'REGION' (first integer is chromosome), 
#   11:'CHR_ID' (also within CHR_POS), 12:'CHR_POS' (chromosome position), 13:'REPORTED GENES', 
#   14:'MAPPED_GENE', 15:'UPSTREAM_GENE_ID', 16:'DOWNSTREAM_GENE_ID', 17:'SNP_GENE_IDS', 
#   18:'UPSTREAM_GENE_DISTANCE', 19:'DOWNSTREAM_GENE_DISTANCE', 20:'STRONGEST SNP-RISK ALLELE' 
#   21:'SNPS', 22:'MERGED', 23:'SNP_ID_CURRENT', 24:'CONTEXT' (variant type), 25:'INTERGENIC', 
#   26:'ALLELE FREQUENCY', 27:'P-VALUE', 28:'PVALUE_MLOG', 29:'P-VALUE (TEXT)', 30:'OR or BETA' 
#   31:'95% CI (TEXT)', 32:'PLATFORM [SNPS PASSING QC]', 
#   33:'CNV' (all N: not a copy number variant) 
gwas_headers = gwas_matrix.pop(0) 

# keep columns of interest - link:5, study:6, trait:7, sample:8, region:10, position:12, genes:13, 
#     mapped genes:14, strongest allele and snp rsid:20, int rsid:23, snp type:24, 
#     allele frequency:26, pvalue:27, log pvalue:28, pvalue text:29, OR BETA:30, 
#     95% confidence interval:31 
gwas_cols = [ 5, 6, 7, 8, 10, 12, 13, 14, 20, 23, 24, 26, 27, 28, 29, 30, 31 ] 
slim_headers = slim_matrix([gwas_headers], gwas_cols) 
slim_headers = slim_headers[0] 
rsid_col = slim_headers.index('SNP_ID_CURRENT') 

slim_gwas = slim_matrix(gwas_matrix, gwas_cols) # new matrix with fewer columns 

# make rsid the first column and sort by rsid 
gwas_rsids = rsid_matrix(slim_gwas, rsid_col) 
gwas_rsids.sort() 

# fix header column with rsid first 
rsid_col = slim_headers.pop(rsid_col) 
slim_headers = [rsid_col] + slim_headers 

# get new header indeces 
allele = slim_headers.index('STRONGEST SNP-RISK ALLELE') 
effect = slim_headers.index('OR or BETA') 
ci = slim_headers.index('95% CI (TEXT)') 

gwas_w_allele = [] 

for row in gwas_rsids: 
    if row[allele][-1] != '?': # remove entries that don't have a specified allele 

        try: # remove entries that don't have a specified beta coefficient 
            effect_size = float(row[effect]) 

            # TODO: move this to a new loop to debug / look at output from just float() line 
            ci_txt = row[ci] 
            if 'crease' not in ci_txt: # odds ratio instead of beta coefficient in OR or beta 
                row[effect] = str( math.log(effect_size) ) # beta coefficient = ln(odds ratio) 
                row[ci] = row[effect] + ' standard deviation unit (odds ratio: ' + ci_txt + ')' 
            else: 
                row[ci] = 'standard deviation ' + ci_txt 
            gwas_w_allele.append(row) 

        except (IndexError, ValueError): # missing the phenotypic effect 
            1 

# add headers back to matrix 
g_len = len(gwas_w_allele) 
gwas_w_allele = [slim_headers] + gwas_w_allele 

# write output matrix to slim_gwas_catalog.txt 
output_file = open('slim_gwas_catalog.txt', 'w', encoding = "utf8") 

print('# Original gwas catalog cut down from 496,340 entries to', str(g_len), 'and ordered by',
      'SNP id\n# This file only contains SNP entries with a specified allele and measured beta', 
      'coefficient\n# 34 original catalog columns have been cut down to 17 columns of interest\n', 
      file = output_file) 

for row in gwas_w_allele: 
    out_s = '' 
    for i in row: 
        if i == '': 
            out_s += 'NR\t' 
        else: 
            out_s += str(i) + '\t' 

    out_s = out_s.strip() 
    out_s += '\n' 
    output_file.write(out_s) 
output_file.close() 

# TODO: use GWAS catalog to create a list of annotated polymorphisms and their rsids 
# Sort catalog by chromosome (REGION), then location (CHR_POS) 

