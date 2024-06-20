# Human_Genome_Analysis
Human_Genome_Analysis uses a GWAS catalog, an annotated database of genetic variants, combined with an individual human’s genetic data, to output the most significant effects for that person’s genetic variants. 

Many Genome-Wide Association Studies exist to characterize genetic variants and how they affect the biological characteristics of a person. These studies are assembled in GWAS catalogs. With the goal of expanding the public knowledge of human genetic diversity, it is important to use a large database of research that is also well-vetted, updating, transparent, and produces the strongest scientific claims possible. This project uses the NIH’s EMBL-EBI GWAS catalog, which contains the single nucleotide variants from every eligible study published on PubMed. For the inclusion criteria, a study must have genome-wide coverage genotyping and analysis of over one hundred thousand SNPs across the genome. The catalog contains one entry, per trait, per publication. The exact file used is ‘all
associations v1.0’ downloaded from https://www.ebi.ac.uk/gwas/docs/file-downloads. The catalog is updated on a weekly basis. Future versions of the catalog may have different structure or inclusion criteria than what is described here. 

An individual person’s “SNP genome” is a data file of single nucleotide variants in that person. For the individual SNP genomes to run as input, I used genomes from those publicly available on snpedia at https://www.snpedia.com/index.php/Genomes. Each row of these files represents an SNP, and contains the SNP’s rsid, chromosome, the position of the nucleotide, and the person’s genotype. GWAS catalog column descriptions: https://www.ebi.ac.uk/gwas/docs/fileheaders.

The code filters the data for variants both in the genome and the GWAS catalog. The SNPs in both files are identified by their rsid, then the annotations of each SNP from the GWAS catalog is saved. The GWAS annotations include the SNP’s phenotype, effect on phenotype, genes, allele, and population frequency. That data processing generates a list of the person’s genetic variants, each with its annotation from the GWAS catalog. These variants are ranked by their significance - the strength / magnitude of the variant’s effect on the phenotype, as based on the beta coefficient. Beta coefficients are reported in units of standard deviations, and beta coefficient = ln(odds ratio). These highest-ranked variants are output as the ‘significant’ variants of the person. 

All input files must be coded in utf-8. Processing also assumes tab-deliminated files. 

## File Structure 

### genomes folder 
Contains SNPome files. 

### file_processing.py 

### gwas_processor.py 

### penotype_predictor.py 

### vcf_genome_reader.py 

## References

### NHGRI-EBI GWAS Catalog
Sollis E, Mosaku A, Abid A, Buniello A, Cerezo M, Gil L, Groza T, Güneş O, Hall P, Hayhurst J, Ibrahim A, Ji Y, John S, Lewis E, MacArthur JAL, McMahon A, Osumi-Sutherland D, Panoutsopoulou K, Pendlington Z, Ramachandran S, Stefancsik R, Stewart J, Whetzel P, Wilson R, Hindorff L, Cunningham F, Lambert SA, Inouye M, Parkinson H, Harris LW. The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource. Nucleic Acids Res. 2022 Nov 9:gkac1010. doi: 10.1093/nar/gkac1010. Epub ahead of print. PMID: 36350656.

### SNPedia Genomes
Genomes. SNPedia. Retrieved May 5, 2023, from https://www.snpedia.com/index.php/Genomes