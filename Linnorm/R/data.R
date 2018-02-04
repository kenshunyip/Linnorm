#' scRNA-seq data from Islam et al. 2011
#' 
#' GEO accession GSE29087:
#' 92 single cells (48 mouse embryonic stem cells, 44 mouse embryonic fibroblasts and 4 negative controls) were analyzed by single-cell tagged reverse transcription (STRT).
#' @docType data
#' @name Islam2011
#' @usage data(Islam2011) 
#' @format A matrix with 22936 rows (genes) and 96 columns (samples). The first 48 columns are ES cells, the following 44 columns are mouse embryonic fibroblasts and the remaining 4 columns and negative controls. Data is in raw counts format.
#' @references Islam et al. (2011) Genome Res 2011 Jul;21(7):1160-7. PMID: 21543516
NULL


#' Partial RNA-seq data from TCGA LIHC (Liver Hepatocellular Carcinoma)
#' 
#' TPM Expression data
#' 
#' @docType data
#' @name LIHC
#' @usage data(LIHC)
#' @format A matrix with 25914 rows (genes) and 10 columns (samples). The first 5 columns are Tumor samples, the remaining 5 columns are adjacent Normal samples. They are paired samples from 5 individuals. Data is in TPM format.
#' @references https://tcga-data.nci.nih.gov/
NULL

#' Partial RNA-seq data from SEQC/MAQC-III Sample A
#' 
#' Raw Count data
#' 
#' @docType data
#' @name SEQC
#' @usage data(SEQC)
#' @format A matrix with 50227 rows (genes) and 10 columns (samples).
#' @references SEQC/MAQC-III Consortium. A comprehensive assessment of RNA-seq accuracy, reproducibility and information content by the Sequencing Quality Control Consortium. Nature biotechnology 32.9 (2014): 903-914.
NULL
