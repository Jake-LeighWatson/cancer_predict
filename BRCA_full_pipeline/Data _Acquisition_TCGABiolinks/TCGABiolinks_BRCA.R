library(SummarizedExperiment)
library(TCGAbiolinks)
library(maftools)
library(dplyr)


##### Tumour Samples ####

#Create query for BRCA  cancer from TCGA, RNA-Seq FPKM values for tumour samples
query_BRCA_tumour <- GDCquery(project = "TCGA-BRCA",
                              data.category = "Transcriptome Profiling", 
                              data.type = "Gene Expression Quantification", 
                              experimental.strategy = "RNA-Seq",
                              workflow.type = "HTSeq - FPKM",
                              sample.type = 'Primary Tumor')

#Download tumour query 
GDCdownload(query_BRCA_tumour)

#Prepare tumour query
BRCA_FPKM_tumour <- GDCprepare(query = query_BRCA_tumour, save = TRUE,save.filename = 'BRCA_tumour.rda')

#Create dataframe containing the tumour FPKM values (stored in assay) for each gene and sample
FPKM_BRCA_tumour <- as.data.frame(SummarizedExperiment::assay(BRCA_FPKM_tumour))

#Create dataframe containing gene symbols (should have same index as the ensembl ids)
gene_symbols_tumour <- data.frame(BRCA_FPKM_tumour@rowRanges@elementMetadata$external_gene_name)
#rename col
colnames(gene_symbols_tumour) <- 'Gene_symbol'

#Join the gene_symbols DF to the FPKM DF
FPKM_BRCA_tumour <- cbind(gene_symbols_tumour,FPKM_BRCA_tumour)


#Write csv for tumour FPKM values
write.csv(FPKM_BRCA_tumour,file = 'FPKM_BRCA_tumour')



#### Normal Samples ####


#Create query for BRCA cancer from TCGA, RNA-Seq FPKM values for normal samples
query_BRCA_normal <- GDCquery(project = "TCGA-BRCA",
                              data.category = "Transcriptome Profiling", 
                              data.type = "Gene Expression Quantification", 
                              experimental.strategy = "RNA-Seq",
                              workflow.type = "HTSeq - FPKM",
                              sample.type = 'Solid Tissue Normal')

#Download normal query
GDCdownload(query_BRCA_normal)

#Prepare normal query
BRCA_FPKM_normal <- GDCprepare(query = query_BRCA_normal, save = TRUE,save.filename = 'BRCA_normal.rda')

#Create dataframe containing the normal FPKM values (stored in assay) for each gene and sample
FPKM_BRCA_normal <- as.data.frame(SummarizedExperiment::assay(BRCA_FPKM_normal))

#Create dataframe containing gene symbols (should have same index as the ensembl ids)
gene_symbols_normal <- data.frame(BRCA_FPKM_normal@rowRanges@elementMetadata$external_gene_name)
#rename col
colnames(gene_symbols_normal) <- 'Gene_symbol'

#Join the gene_symbols DF to the FPKM DF
FPKM_BRCA_normal <- cbind(gene_symbols_normal,FPKM_BRCA_normal)


#write normal csv
write.csv(FPKM_BRCA_normal,file = 'FPKM_BRCA_normal')

