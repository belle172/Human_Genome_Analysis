# -*- coding: utf-8 -*- 
''' 
Python 3.10 
Created Mar 30 2023 
Jasper Bellefeuille - jasperbellefeuille@gmail.com 
Repository: Human_Genome_Analysis/phenotype_predictor.py 

This script takes an individual's genotype file and outputs the strongest phenotypic effects 

Assumptions: 
    slim_gwas_catalog.txt has been created by first running gwas_processor.py 
    The SNP genome files are in the path 
    An individual's genome file contains rsids, or location based on hg38 as reference genome 
''' 

from matplotlib import pyplot 
from file_processing import get_tsv_matrix, rsid_matrix 

def genome_processing(genome_matrix, name): 

    # Turn genome file into matrix where rsid is the first column 
    genome_rsids = rsid_matrix(genome_matrix, 0) 
    genome_rsids.sort() # sort the genome matrix by rsid, so it matches the GWAS catalog 
    geno_len = len(genome_rsids) # number of genetic locations genotyped in the genome 

    # Initialize processed gwas catalog 
    gwas_rsids = get_tsv_matrix('slim_gwas_catalog.txt') 
    gwas_hdrs = gwas_rsids.pop(0) # headers from first line of catalog 
    gwas_rsids = rsid_matrix(gwas_rsids, 0) # reorder catalog so rsid is also the first column 
    gwas_len = len(gwas_rsids) 

    # rsid_num chromosome genome_position genotype gwas_headers 
    gnome_hdrs = [gwas_hdrs[0]] + ['chromosome', 'genome position', 'genotype'] + gwas_hdrs[1:] 

    annotated_snps = [] # List of gwas entries for the annotated SNPs 

    # format: rsid, chromosome, position, genotype, rest of gwas columns 
    for snp in genome_rsids: # for loop takes about 50 seconds to run 

        # current genome snp is higher than the next cataloged snp, so go to next catalog snp 
        while snp[0] > gwas_rsids[0][0]: 
            gwas_rsids.pop(0) 

        # After possibly moving up in the catalog, check if the next snp in the genome is in the 
        while snp[0] == gwas_rsids[0][0]: # gwas catalog (might be multiple entries of SNP) 
            gwas_line = gwas_rsids.pop(0) # go to the next snp in the catalog 
            gwas_line = snp + gwas_line[1:] 
            annotated_snps.append(gwas_line) # append the matching snp to the list 
        # else: current genome snp < next cataloged snp, meaning its not cataloged 

    # get indeces of columns in the GWAS catalog matrix 
    allele = gnome_hdrs.index('STRONGEST SNP-RISK ALLELE') 
    geno   = gnome_hdrs.index('genotype') 
    effect = gnome_hdrs.index('OR or BETA') 

    measured_snps = [] 

    # Get just the snps that have an allele in the person's genotype 
    for snp in annotated_snps: 
        if (snp[allele][-1]) in snp[geno]: 
            effect_size = float(snp[effect]) 
            line = [ effect_size ] + snp[0:effect] + snp[ (effect+1): ] 
            measured_snps.append(line) 

    measured_snps.sort(reverse = True) # sort SNPs by the effect size / beta coefficient 

    effect = gnome_hdrs.pop(effect) 
    gnome_hdrs = [effect] + gnome_hdrs 
    effect_output = [gnome_hdrs] + measured_snps 

    # write the snps ordered by effect size 
    f_name = name + '_snps_by_effect_size.txt' 
    output_file = open(f_name, 'w', encoding = 'utf8') 
    print('# File contains', str(len(measured_snps)), 'SNP match entries in the GWAS catalog', 
          'that have a measured phenotypic effect size, ordered by largest effect size first\n#', 
          'Original catalog contains', str(gwas_len), 'entries, and the genome contains', 
          str(geno_len), ' SNPs\n', file = output_file) 
    
    for row in effect_output: 
        out_s = '' 
        for i in row: 
            out_s += str(i) + '\t' 
        out_s = out_s.strip() 
        out_s += '\n' 
        output_file.write(out_s) 
    output_file.close() 

    return effect_output 


# ============================================================================= 
# run get_tsv_matrix() and genome_processing() on input genome 
# public snp genomes from snpedia - https://www.snpedia.com/index.php/Genomes 
# ============================================================================= 

# TODO: update so any SNPome file can be input instead of manually naming it in the two lines 
# genome files are organized by chromosome, then position in chromosome. columns - rsid, 
lilly = get_tsv_matrix('genome_Lilly_Mendel_v4.txt') # chromosome, position, genotype 
measured_snps = genome_processing(lilly, 'lilly') # function takes a minute 
# denisova = get_tsv_matrix('denisova.23andme.112') 
# measured_snps = genome_processing(denisova, 'denisova') 

genome_hdrs = measured_snps.pop(0) 

# Input percent of the most significant SNPs to display 
percent = 0.001 
cutoff = range( int( percent * len(measured_snps) ) ) 

# =============================================================================
# Write significant data file and phenotype graph 
# ============================================================================= 
s = 'outputting data for the highest ' + str(percent*100) + '% of SNPs by phenotype effect size ' 
s += 'in this genome\nthis list contains ' + str(len(cutoff)) + ' SNPs\n\n' 

out_f = open('significant_snps.txt', 'w', encoding = "utf8") 
out_f.write(s) 

for i in cutoff: 
    s = 'SNP id rs' + str(measured_snps[i][1]) + ', rank ' + str(i+1) + ' in effect size\n' 

    # print the genes associated with the snp 
    if (measured_snps[i][11] == measured_snps[i][12]) or (measured_snps[i][11] == 'NR'): 
        g = measured_snps[i][12] 
    elif measured_snps[i][12] == 'NR': 
        g = measured_snps[i][11] 
    else: 
        g = measured_snps[i][11] + ', ' + measured_snps[i][12] 

    if (',' in g) or ('-' in g): 
        g = 'genes: ' + g 
    else: 
        g = 'gene: ' + g 
    s += g + '\n' 

    s += 'associated phenotype: ' + measured_snps[i][7] 
    ci = measured_snps[i][19] # output the 95% confidence interval if there is one 
    if ci[:9] == '[NR] unit': 
        ci = ci[10:] 
    elif ci[:4] == '[NR]': 
        ci = ci[5:] 

    if measured_snps[i][18][1:-1] == 'AA': 
        s += ', with ' + ci + ' in ' + measured_snps[i][18][1:-1] 
    elif measured_snps[i][18] != 'NR': # if the phenotype effect is specified 
        s += ', with ' + ci + ' in ' + measured_snps[i][18][1:-1] 
    else: 
        s += ' ' + ci 
    s += ' compared to phenotype of reference genome\n' 

    s += 'variant allele associated with phenotype: ' + measured_snps[i][13][-1] + ', population ' 
    s += 'frequency: ' + measured_snps[i][15] + '\ngenotype of inputted genome: ' 
    s += measured_snps[i][4] + '\n' 
    
    s += 'beta coefficient - strength of effect on phenotype: ' 
    s += str(measured_snps[i][0]) + '\n' 

    s += 'chromosome region: ' + measured_snps[i][9] + '\nreference position: ' 
    s += measured_snps[i][10] + ', position on input genome: ' + measured_snps[i][3] + '\n' 

    s += '\n' 
    out_f.write(s) 
out_f.close() 

# output total effect on phenotypes 
pheno = genome_hdrs.index('DISEASE/TRAIT') 
direction = genome_hdrs.index('95% CI (TEXT)') 

traits = {} 
for snp in measured_snps: 
    if 'increase' in snp[direction]: 
        effect = snp[0] 
    elif 'decrease' in snp[direction]: 
        effect = -1 * snp[0] 

    if snp[pheno] not in traits: 
        traits[ snp[pheno] ] = effect 
    else: 
        traits[ snp[pheno] ] += effect 

phenotypes = [] 
total = 0 

# calculate whether each trait was in total decreased or increased 
for (trait, effect) in traits.items(): 
    if effect < 0: 
        crease = 'decrease in ' + trait 
        phenotypes.append( [-1 * effect, crease] ) 
        total += (-1 * effect) 
    else: 
        crease = 'increase in ' + trait 
        phenotypes.append( [effect, crease] ) 
        total += effect 

phenotypes.sort(reverse = True) 

effects, phenos = [], [] 
portion, other = 0, 0 
for (effect, trait) in phenotypes: 
    if portion < (total * 0.7): 
        effects.append(effect/total) 
        phenos.append(trait) 
        portion += effect 
    else: 
        other += (effect/total) 

effects.append(other) 
phenos.append('other') 

def pie_labels(value, value_list): 
    return str( round(value, 1) ) + '%' 

figure, axes = pyplot.subplots(figsize = (12, 8)) 
axes.pie( effects, autopct = lambda value: pie_labels(value, effects), 
                             labels = phenos, startangle = 320, 
                             textprops = dict(color = 'black') ) 
axes.set_title('Portion of total effect for traits in input genome', fontsize = 30) 

