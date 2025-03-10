# Load in catalog 
View(gwas_catalog_v1.associations_e109) 

# Get number of unique studies 
length(unique(gwas_catalog_v1.associations_e109$PUBMEDID)) 