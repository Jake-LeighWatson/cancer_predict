setwd("~/OneDrive - University of Glasgow/Project/THCA/R_data")

load('THCA_FPKM.rdata')

library(SummarizedExperiment)
library(TCGAbiolinks)
library(maftools)
library(dplyr)


##### Tumour Samples ####

#Create query for THCA  cancer from TCGA, RNA-Seq FPKM values for tumour samples
query_THCA_tumour <- GDCquery(project = "TCGA-THCA",
                              data.category = "Transcriptome Profiling", 
                              data.type = "Gene Expression Quantification", 
                              experimental.strategy = "RNA-Seq",
                              workflow.type = "HTSeq - FPKM",
                              sample.type = 'Primary Tumor')

#Download tumour query 
GDCdownload(query_THCA_tumour)

#Prepare tumour query
THCA_FPKM_tumour <- GDCprepare(query = query_THCA_tumour, save = TRUE,save.filename = 'THCA_tumour.rda')

#Create dataframe containing the tumour FPKM values (stored in assay) for each gene and sample
FPKM_THCA_tumour <- as.data.frame(SummarizedExperiment::assay(THCA_FPKM_tumour))

#Create dataframe containing gene symbols (should have same index as the ensembl ids)
gene_symbols_tumour <- data.frame(THCA_FPKM_tumour@rowRanges@elementMetadata$external_gene_name)
#rename col
colnames(gene_symbols_tumour) <- 'Gene_symbol'

#Join the gene_symbols DF to the FPKM DF
FPKM_THCA_tumour <- cbind(gene_symbols_tumour,FPKM_THCA_tumour)


#Write csv for tumour FPKM values
write.csv(FPKM_THCA_tumour,file = 'FPKM_THCA_tumour')



#### Normal Samples ####


#Create query for THCA cancer from TCGA, RNA-Seq FPKM values for normal samples
query_THCA_normal <- GDCquery(project = "TCGA-THCA",
                              data.category = "Transcriptome Profiling", 
                              data.type = "Gene Expression Quantification", 
                              experimental.strategy = "RNA-Seq",
                              workflow.type = "HTSeq - FPKM",
                              sample.type = 'Solid Tissue Normal')

#Download normal query
GDCdownload(query_THCA_normal)

#Prepare normal query
THCA_FPKM_normal <- GDCprepare(query = query_THCA_normal, save = TRUE,save.filename = 'THCA_normal.rda')

#Create dataframe containing the normal FPKM values (stored in assay) for each gene and sample
FPKM_THCA_normal <- as.data.frame(SummarizedExperiment::assay(THCA_FPKM_normal))

#Create dataframe containing gene symbols (should have same index as the ensembl ids)
gene_symbols_normal <- data.frame(THCA_FPKM_normal@rowRanges@elementMetadata$external_gene_name)
#rename col
colnames(gene_symbols_normal) <- 'Gene_symbol'

#Join the gene_symbols DF to the FPKM DF
FPKM_THCA_normal <- cbind(gene_symbols_normal,FPKM_THCA_normal)


#write normal csv
write.csv(FPKM_THCA_normal,file = 'FPKM_THCA_normal')



#### MAF file ####

#Download maf file using GDCquery_Maf
maf_THCA <- GDCquery_Maf("THCA", pipelines = "muse")

#Write the maf file to csv
write.csv(maf_THCA,file = 'MAF_THCA')

#Visualize MAF file data
getSampleSummary(x=read.maf(maf = maf_THCA))
plotmafSummary(maf = read.maf(maf = maf_THCA), rmOutlier = TRUE, addStat = 'mean', dashboard = TRUE)

help(plotmafSummary)

save.image('THCA_FPKM.rdata')
