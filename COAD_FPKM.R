setwd("~/OneDrive - University of Glasgow/Project - Practice/TCGA_R_FPKM")


library(SummarizedExperiment)
library(TCGAbiolinks)


##### Tumour Samples ####

#Create query for COAD cancer from TCGA, RNA-Seq FPKM values for tumour samples
query_COAD_tumour <- GDCquery(project = "TCGA-COAD",
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           experimental.strategy = "RNA-Seq",
                           workflow.type = "HTSeq - FPKM",
                           sample.type = 'Primary Tumor')

#Download tumour query 
GDCdownload(query_COAD_tumour)

#Prepare tumour query
coad_FPKM_tumour <- GDCprepare(query = query_COAD_tumour, save = TRUE,save.filename = 'coad_tumour.rda')
coad_FPKM_tumour

#Create dataframe containing the tumour FPKM values (stored in assay) for each gene and sample
FPKM_COAD_tumour <- as.data.frame(SummarizedExperiment::assay(coad_FPKM_tumour))
FPKM_COAD

#Write csv for tumour FPKM values
write.csv(FPKM_COAD_tumour,file = 'FPKM_COA_tumour')



#### Normal Samples ####


#Create query for COAD cancer from TCGA, RNA-Seq FPKM values for normal samples
query_COAD_normal <- GDCquery(project = "TCGA-COAD",
                              data.category = "Transcriptome Profiling", 
                              data.type = "Gene Expression Quantification", 
                              experimental.strategy = "RNA-Seq",
                              workflow.type = "HTSeq - FPKM",
                              sample.type = 'Solid Tissue Normal')

#Download normal query
GDCdownload(query_COAD_normal)

#Prepare normal query
coad_FPKM_normal <- GDCprepare(query = query_COAD_normal, save = TRUE,save.filename = 'coad_normal.rda')
coad_FPKM_normal

#Create dataframe containing the normal FPKM values (stored in assay) for each gene and sample
FPKM_COAD_normal <- as.data.frame(SummarizedExperiment::assay(coad_FPKM_normal))

#write normal csv
write.csv(FPKM_COAD_normal,file = 'FPKM_COAD_normal')


