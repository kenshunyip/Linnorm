#Linnorm

Linnorm is an R package for the analysis of RNA-seq, scRNA-seq, ChIP-seq count data or any large scale count data. It transforms such datasets for parametric tests. In addition to the transformtion function (Linnorm), the following pipelines are implemented: 1. Library size/Batch effect normalization (Linnorm.Norm), 2. Cell subpopluation analysis and visualization using t-SNE or PCA K-means clustering or Hierarchical clustering (Linnorm.tSNE, Linnorm.PCA, Linnorm.HClust), 3. Differential expression analysis or differential peak detection using limma (Linnorm.limma), 4. Highly variable gene discovery and visualization (Linnorm.HVar), 5. Gene correlation network analysis and visualization (Linnorm.Cor), 6. Stable gene selection for scRNA-seq data; for users without or do not want to rely on spike-in genes (Linnorm.SGenes). 7. Data imputation. (under development) (Linnorm.DataImput). Linnorm can work with raw count, CPM, RPKM, FPKM and TPM. Additionally, the RnaXSim function is included for simulating RNA-seq data for the evaluation of DEG analysis methods.

#Installation:


source("https://bioconductor.org/biocLite.R")

biocLite("Linnorm")


#Manual:


Please see Linnorm_User_Manual.pdf