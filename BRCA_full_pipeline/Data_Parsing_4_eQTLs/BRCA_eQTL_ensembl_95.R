setwd("~/OneDrive - University of Glasgow/Project/BRCA/Machine_learning_final/Pipeline/Data_parsing_4_eQTLs")

library(biomaRt)
library(dplyr)
library(plyr)

#Load eQTL data for cancer type
cis_eQTL_95 <- read.csv('BRCA_cis_eQTL_95.csv',sep = ',')
trans_eQTL_95 <- read.csv('BRCA_trans_eQTL_95.csv',sep = ',')

#Remove added X col
cis_eQTL_95 <- select(cis_eQTL_95,-X)
trans_eQTL_95 <- select(trans_eQTL_95,-X)

#Get gene entrez IDs from eQTL tables
cis_95_genes_entrez <- cis_eQTL_95$Entrez_ID
trans_95_genes_entrez <- trans_eQTL_95$Entrez_ID

#Use human ensembl 
ensemblHuman = useDataset("hsapiens_gene_ensembl",mart=useMart("ensembl"))

#getBM via entrez
cis_95_entrez <- getBM(filters= "entrezgene_id", attributes= c("hgnc_symbol","ensembl_gene_id",'entrezgene_id',"chromosome_name"),values=cis_95_genes_entrez,mart= ensemblHuman)
trans_95_entrez <- getBM(filters= "entrezgene_id", attributes= c("hgnc_symbol","ensembl_gene_id",'entrezgene_id',"chromosome_name"),values=trans_95_genes_entrez,mart= ensemblHuman)

#Drop rows which contain genes not on a main chromosome => drop duplicated genes
cis_95_entrez <- subset(cis_95_entrez, nchar(as.character(cis_95_entrez$chromosome_name)) <= 2)
trans_95_entrez <- subset(trans_95_entrez, nchar(as.character(trans_95_entrez$chromosome_name)) <= 2)


#Repeat for gene symbol 

cis_95_genes_symbol <- cis_eQTL_95$Gene_symbol
trans_95_genes_symbol <- trans_eQTL_95$Gene_symbol

cis_95_symbol <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","ensembl_gene_id",'entrezgene_id',"chromosome_name"),values=cis_95_genes_symbol,mart= ensemblHuman)
trans_95_symbol<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","ensembl_gene_id",'entrezgene_id',"chromosome_name"),values=trans_95_genes_symbol,mart= ensemblHuman)

#Drop rows which contain genes not on a main chromosome => drop duplicated genes
cis_95_symbol <- subset(cis_95_symbol, nchar(as.character(cis_95_symbol$chromosome_name)) <= 2)
trans_95_symbol <- subset(trans_95_symbol, nchar(as.character(trans_95_symbol$chromosome_name)) <= 2)


#Create table joining original cis_eQTL table with converted tables (entrez and symbol)
cis_joined_entrez <- left_join(cis_eQTL_95, cis_95_entrez, by = c("Entrez_ID"="entrezgene_id"))
cis_joined_symbol <- left_join(cis_eQTL_95, cis_95_symbol, by = c("Gene_symbol"="hgnc_symbol"))

trans_joined_entrez <- left_join(trans_eQTL_95, trans_95_entrez, by = c("Entrez_ID"="entrezgene_id"))
trans_joined_symbol <- left_join(trans_eQTL_95, trans_95_symbol, by = c("Gene_symbol"="hgnc_symbol"))

#Write tables to csv for 95 cis and trans
write.table(cis_joined_entrez, file = 'BRCA_cis_entrez_95.csv', sep = '\t')
write.table(cis_joined_symbol, file = 'BRCA_cis_symbol_95.csv', sep = '\t')

write.table(trans_joined_entrez, file = 'BRCA_trans_entrez_95.csv', sep = '\t')
write.table(trans_joined_symbol, file = 'BRCA_trans_symbol_95.csv', sep = '\t')


