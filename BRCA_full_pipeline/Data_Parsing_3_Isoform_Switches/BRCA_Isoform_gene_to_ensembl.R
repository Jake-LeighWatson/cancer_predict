setwd("~/OneDrive - University of Glasgow/Project/BRCA/Machine_learning_final/Pipeline/Data_parsing_3_iso")

library(biomaRt)
library(dplyr)

#Load isoform_switch data for cancer type
isoform_swiches <- read.csv('BRCA_isoform_switches.csv',sep = ',')

#Remove added X col
isoform_swiches <- select(isoform_swiches,-X)

#Get genes from isoform_switch table
genes <- isoform_swiches$Entrez_ID

#Use human ensembl - GRCh38
ensemblHuman = useDataset("hsapiens_gene_ensembl",mart=useMart("ensembl"))
isoform_switches_ensembl_symbol <- getBM(filters= "entrezgene_id", attributes= c("hgnc_symbol","ensembl_gene_id",'entrezgene_id',"chromosome_name"),values=genes,mart= ensemblHuman)

#Drop rows which contain genes not on a main chromosome => drop duplicated genes
isoform_switches_ensembl_symbol <- subset(isoform_switches_ensembl_symbol, nchar(as.character(isoform_switches_ensembl_symbol$chromosome_name)) <= 2)

#Change hgnc col name to 'Gene_symbol'
colnames(isoform_switches_ensembl_symbol)[1] <- 'Gene_symbol'

write.table(isoform_switches_ensembl_symbol,file = 'BRCA_isoform_switches_ensembl.csv',sep = '\t')

  