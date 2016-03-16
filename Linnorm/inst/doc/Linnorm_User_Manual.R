## ---- echo=FALSE, eval=TRUE----------------------------------------------
library("knitr")
opts_chunk$set(tidy=FALSE,results="asis", message=FALSE,warning=FALSE)


## ---- echo=TRUE----------------------------------------------------------
library(seqc)

SampleA <- ILM_aceview_gene_BGI[,grepl("A_",colnames(ILM_aceview_gene_BGI))]
rownames(SampleA) <- ILM_aceview_gene_BGI[,2]

SampleB <- ILM_aceview_gene_BGI[,grepl("B_",colnames(ILM_aceview_gene_BGI))]
rownames(SampleB) <- ILM_aceview_gene_BGI[,2]

## ---- echo=TRUE----------------------------------------------------------
set.seed(12345)
SampleA3 <- SampleA[,sample(1:80,6)]
SampleB3 <- SampleB[,sample(1:80,8)]


## ---- echo=TRUE----------------------------------------------------------
datamatrix <- cbind(SampleA3,SampleB3)

## ---- echo=TRUE----------------------------------------------------------
design <- c()
for (i in 1:6) {
	design <- c(design, 1)
}
for (i in 1:8) {
	design <- c(design, 2)
}
design <- model.matrix(~ 0+factor(design))
colnames(design) <- c("group1", "group2")
rownames(design) <- colnames(datamatrix)

## ---- echo=TRUE----------------------------------------------------------
library(Linnorm)
DEG_Results <- Linnorm.limma(datamatrix,design)
#The DEG_Results matrix contains DEG analysis results.

## ---- echo=TRUE----------------------------------------------------------
library(Linnorm)
BothResults <- Linnorm.limma(datamatrix,design,output="Both")

#To separate results into two matrices for analysis:
DEG_Results <- BothResults[[1]]
#Or
DEG_Results <- BothResults$DEResults
#The DEG_Results matrix now contains DEG analysis results.

TransformedMatrix <- BothResults[[2]]
#Or
TransformedMatrix <- BothResults$TMatrix
#The TransformedMatrix matrix now contains a Linnorm Normalized dataset.

## ---- eval=FALSE---------------------------------------------------------
#  write.table(DEG_Results, "DEG_Results.txt", quote=F, sep="\t", col.names=TRUE,
#  row.names=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  Genes15 <- DEG_Results[order(DEG_Results[,"adj.P.Val"]),][1:15,]
#  print(Genes15)

## ----kable, echo=FALSE---------------------------------------------------
Genes15 <- DEG_Results[order(DEG_Results[,"adj.P.Val"]),][1:15,]
library(knitr)
kable(Genes15, digits=4)

## ---- echo=TRUE, fig.height=6, fig.width=6-------------------------------
#Remove Genes which fold changes are INF. You can skip this if there is no INF 
#values in the fold change column.
NoINF <- DEG_Results[which(!is.infinite(DEG_Results[,1])),]

#Record significant genes for coloring
SignificantGenes <- NoINF[NoINF[,5] <= 0.05,1]

#Draw volcano plot with Log q values. 
#Green dots are non-significant, red dots are significant.
plot(x=NoINF[,1], y=NoINF[,5], col=ifelse(NoINF[,1] %in% SignificantGenes, 
"red", "green"),log="y", ylim = rev(range(NoINF[,5])), main="Volcano Plot", 
xlab="log2 Fold Change", ylab="q values", cex = 0.5)

## ---- echo=TRUE----------------------------------------------------------
library(seqc)
SampleA <- ILM_aceview_gene_BGI[,grepl("A_",colnames(ILM_aceview_gene_BGI))]
rownames(SampleA) <- ILM_aceview_gene_BGI[,2]



## ---- echo=TRUE----------------------------------------------------------
library(Linnorm)
#We are only normalizing part of the matrix for an example.
SampleA <- SampleA[,1:20]
Normalized <- Linnorm(SampleA)

## ---- eval=FALSE---------------------------------------------------------
#  #You can use this file with Excel.
#  write.table(Normalized, "Normalized_Matrix.txt", quote=F, sep="\t",
#  col.names=TRUE, row.names=TRUE)

## ---- echo=TRUE----------------------------------------------------------
Normalized <- Normalized[rowSums(Normalized == 0) <= length(Normalized[1,])/4,]

## ---- echo=TRUE----------------------------------------------------------
#Here, we randomly choose 100 genes for an example. 
#However, users can provide a custom gene list here.
set.seed(12345)
Normalized <- Normalized[sample(1:length(Normalized[,1]),100),]

#Calculate correlation.
PCC <- cor(t(Normalized))

## ---- echo=TRUE, fig.height=8, fig.width=8-------------------------------
#Please install these libraries if you need to.
library(RColorBrewer)
library(gplots)
heatmap.2(as.matrix(PCC), Rowv=TRUE, Colv=TRUE, dendrogram='none',
symbreaks=TRUE , trace="none", xlab = 'Genes', ylab = "Genes", 
density.info='none' ,key.ylab=NA, col = colorRampPalette(c("blue",  "white", 
"yellow"))(n = 1000), lmat=rbind(c(4, 3), c(2, 
1)),cexRow=0.3,cexCol=0.3,margins = c(8, 8))

## ---- echo=TRUE----------------------------------------------------------
library(seqc)

SampleA <- ILM_aceview_gene_BGI[,grepl("A_",colnames(ILM_aceview_gene_BGI))]
rownames(SampleA) <- ILM_aceview_gene_BGI[,2]

SampleB <- ILM_aceview_gene_BGI[,grepl("B_",colnames(ILM_aceview_gene_BGI))]
rownames(SampleB) <- ILM_aceview_gene_BGI[,2]

## ---- echo=TRUE----------------------------------------------------------
datamatrix <- cbind(SampleA,SampleB)

## ---- echo=TRUE----------------------------------------------------------
library(Linnorm)
LNormData <- Linnorm(datamatrix)

## ---- echo=TRUE----------------------------------------------------------
#Let LNormData be a matrix of Linnorm Normalized dataset. 
Newdata <- exp(LNormData)

## ---- echo=TRUE----------------------------------------------------------
#Now, we can calculate fold changes between sample set 1 and sample set 2.
#Index of sample set 1 from LNormData and Newdata:
set1 <- 1:80
#Index of sample set 2 from LNormData and Newdata:
set2 <- 81:160

#Define a function that calculates log 2 fold change:
log2fc <- function(x) {
return(log(mean(x[set1])/mean(x[set2]),2))
}

#Calculate log 2 fold change of each gene in the dataset:
foldchanges <- unlist(apply(Newdata,1,log2fc))

#Put resulting data into a matrix
FCMatrix <- matrix(nrow=length(foldchanges),ncol=1)
rownames(FCMatrix) <- rownames(LNormData)
colnames(FCMatrix) <- c("Log 2 Fold Change")
FCMatrix[,1] <- foldchanges

#Now FCMatrix contains fold change results.

## ---- echo=TRUE, fig.height=6, fig.width=6-------------------------------
Density <- density(foldchanges)
plot(Density$x,Density$y,type="n",xlab="Log 2 Fold Change", ylab="Probability Density",)
lines(Density$x,Density$y, lwd=1.5, col="blue")
title("Probability Density of Fold Change.\nSEQC Sample A vs Sample B")
legend("topright",legend=paste("mean = ", round(mean(foldchanges),2),
"\nStdev = ", round(sd(foldchanges),2)))

## ---- echo=TRUE----------------------------------------------------------
library(seqc)

SampleA <- ILM_aceview_gene_BGI[,grepl("A_",colnames(ILM_aceview_gene_BGI))]
rownames(SampleA) <- ILM_aceview_gene_BGI[,2]

## ---- echo=TRUE----------------------------------------------------------
library(Linnorm)
#This will generate two sets of RNA-seq data with 3 replicates each. 
#It will have 20000 genes totally with 5000 genes being differentially 
#expressed. It has the Poisson distribution.
SimulatedData <- RnaXSim(SampleA)

## ---- echo=TRUE----------------------------------------------------------
#Index of differentially expressed genes.
DE_Index <- SimulatedData[[2]]

#Expression Matrix
ExpMatrix <- SimulatedData[[1]]

#Sample Set 1
Sample1 <- ExpMatrix[,1:3]

#Sample Set 2
Sample2 <- ExpMatrix[,4:6]

## ---- echo=TRUE----------------------------------------------------------
library(seqc)

SampleA <- ILM_aceview_gene_BGI[,grepl("A_",colnames(ILM_aceview_gene_BGI))]
rownames(SampleA) <- ILM_aceview_gene_BGI[,2]

## ---- echo=TRUE----------------------------------------------------------
#Number of replicates for each sample set.
NumReplicates <- 5
#Number of Genes in the Samples.
NumGenes <- 5000
#Number of Differentially Expressed Genes.
NumDiffGenes <- 1000
#Distribution. Put "NB" for Negative Binomial, "Gamma" for Gamma,
#"Poisson" for Poisson and "LogNorm" for Log Normal distribution.
Distribution <- "Gamma"

## ---- echo=TRUE----------------------------------------------------------
library(Linnorm)
SimulatedData <- RnaXSim(SampleA, distribution=Distribution,
NumRep=NumReplicates, NumDiff = NumDiffGenes, NumFea = NumGenes)

## ---- echo=TRUE----------------------------------------------------------
#Index of differentially expressed genes.
DE_Index <- SimulatedData[[2]]

#Expression Matrix
ExpMatrix <- SimulatedData[[1]]

#Sample Set 1
Sample1 <- ExpMatrix[,1:3]

#Sample Set 2
Sample2 <- ExpMatrix[,4:6]

## ---- echo=TRUE----------------------------------------------------------
#Obtain a normalized dataset for an example.
library(seqc)
SampleA <- ILM_aceview_gene_BGI[,grepl("A_",colnames(ILM_aceview_gene_BGI))]
rownames(SampleA) <- ILM_aceview_gene_BGI[,2]
SampleA <- SampleA[,1:4]
library(Linnorm)
LNormData <- Linnorm(SampleA)

#LNormData is a matrix of Linnorm Normalized dataset.
#XPMdata is a CPM or TPM matrix.
XPMdata <- ((exp(LNormData) - 1)/mean(colSums(LNormData))) * 1000000
#Now, XPMdata contains a CPM dataset if the original data is raw count or CPM.
#It contains a TPM dataset if the original data is RPKM, FPKM or TPM.


#On the other hand, if LNormData was filtered or changed in some unknown ways;
# and that there is no way to obtain the original dataset for the calculation
#of CPM or TPM. Then, you can convert it this way:
XPMdata <- exp(LNormData) - 1
for (i in seq_along(XPMdata[1,])) {
	XPMdata[,i] <- (XPMdata[,i]/sum(XPMdata[,i])) * 1000000
}
#Now, XPMdata contains a CPM dataset if the original data is raw count or CPM. 
#It contains a TPM dataset if the original data is RPKM, FPKM or TPM.


