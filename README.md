## MSc in Precision Medicine 2021

###  Cancer types stratification based on genomics and transcriptomics data (i.e. genomic variants, gene expression and isoform switching).

### Technology: Bioinformatics for precision medicine
Machine learning classifiers based on genomic variants (SNPs variants within eQTL regions; protein-coding variants), isoform switching, gene expression data to identify important features discriminating between distinct cancer types.

### Databases: The Cancer Genome Atlas (TCGA), GTEx, PancanQTL, iso-kTSP.
* TCGA: Cancer samples from 33 cancer types. https://portal.gdc.cancer.gov/
* GTEx: Gene expression (TPM) and protein-coding variants (per gene) from GTEx.
https://gtexportal.org/
* PancanQTL: Cis-/trans-eQTLs in TCGA samples (regulatory variants). 
http://bioinfo.life.hust.edu.cn/PancanQTL/
http://gong_lab.hzau.edu.cn/PancanQTL/download
Access to all PancanQTL data from https://www.synapse.org/#!Synapse:syn11305829/files/
* iso-kTSP: Isoform switches in TCGA samples.
https://pubmed.ncbi.nlm.nih.gov/25578962/

### Input data: extraction of features from existing databases
* PancanQTL database: Generate matrices (M1, M2) summarising numbers of local (cis-eQTL; M1) and distant (trans-eQTLs; M2) expression quantitative trait loci (eQTL) per gene.
* GTEx/TCGA database: Generate matrix (M3) summarising numbers of protein-coding variants per gene. Extract matrix (referred to as M4) summarising TPM values per gene per sample.
* iso-kTSP resource: Extract genes associated with isoform switches and refine feature space based on that gene subset. Generate a subset matrix (M5) from M4 selecting only those genes associated with known isoform switches. Another option would be to also include the single isoform genes as part of a subset matrix (M6), as those genes may well drive cancer type discrimination, but they would never be associated with isoform switches by definition (due to the presence of only one isoform/gene – no option to switch).

### Machine learning models & cancer types prediction
Build machine learning classifiers (e.g. Random Forest) based on each matrix and compare model performances in classifying cancer types. This step will include dividing datasets into training/test sets, performing cross-validation and evaluating performance metrics on unseen data. Additional work would include multiple input data types (multiple matrices) to build models based on multi-omics data.
