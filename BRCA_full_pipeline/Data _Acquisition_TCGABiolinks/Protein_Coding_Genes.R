
# -- Install Biomart library
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("biomaRt")

# --Reference & manual:
# https://bioconductor.org/packages/release/bioc/html/biomaRt.html
# https://www.bioconductor.org/packages/devel/bioc/manuals/biomaRt/man/biomaRt.pdf

library(biomaRt)

print("Querying Biomart for protein coding genes")
ensemblHuman = useDataset("hsapiens_gene_ensembl",mart=useMart("ensembl"))
humanProteinCodingGenes = getBM(attributes=c("ensembl_gene_id", "external_gene_name", "description"), 
                                filters='biotype', values=c('protein_coding'), mart=ensemblHuman)

print("Save protein coding genes to file")
outputFile = 'ensembl_GRCH38_protein_coding_genes.tsv'
write.table(humanProteinCodingGenes, file=outputFile, quote=FALSE, sep='\t', col.names = NA)
