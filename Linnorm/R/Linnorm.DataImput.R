#' Linnorm Data Imputation Function. (In development)
#'
#' This function performs data imputation for (sc)RNA-seq expression data or large scale count data. It will treat every zero count in the dataset as missing data and replace them with predicted values.
#' @param datamatrix	The matrix or data frame that contains your dataset. It is only compatible with log transformed datasets.
#' @param RowSamples	Logical. In the datamatrix, if each row is a sample and each row is a feature, set this to TRUE so that you don't need to transpose it. Defaults to FALSE.
#' @param showinfo	Logical. Show algorithm running information. Defaults to FALSE.
#' @param MZP	Double >=0, <= 1. Minimum non-Zero Portion Threshold for this function. Genes not satisfying this threshold will be removed. For exmaple, if set to 0.3, genes without at least 30 percent of the samples being non-zero will be removed. Defaults to 0.25.
#' @param LC_F	Double >= 0.01, <= 0.95 or Character "Auto". Filter this portion of the lowest expressing genes. It can be determined automatically by setting to "Auto". Defaults to "Auto".
#' @param max_LC_F	Double >=0, <= 0.95. When LC_F is set to auto, this is the maximum threshold that Linnorm would assign. Defaults to 0.75.
#' @param FG_Recov	Double >=0, <= 1. In the low count gene filtering algorithm, recover this portion of genes that are filtered. Defaults to 0.5.
#' @param method	Character. Method for calculating the distance matrix. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "correlation", "spearman" or "kendall". Any unambiguous substring can be given. Defaults to "euclidean".
#' @param VarPortion	Double >0, <=0.95. Portion of the variance from PCA to be used for data imputation. Defaults to 0.5.
#' @param ... place holder for any new arguments.
#' @details  This function performs data imputation on the dataset. It first generates a distance matrix using principal components from PCA. Then, by default, using the distance matrix as weight, it predicts missing values from each gene using inverse euclidean distance weighted mean.
#' @return This function returns a data matrix.
#' @keywords Linnorm single cell RNA-seq data imputation missing value
#' @export
#' @examples
#' #Obtain example matrix:
#' data(Islam2011)
#' #Transformation:
#' Transformed <- Linnorm(Islam2011)
#' #Data imputation
#' DataImput <- Linnorm.DataImput(Transformed)
#' @import
#' Rcpp
#' RcppArmadillo

Linnorm.DataImput <- function(datamatrix, RowSamples = FALSE, showinfo = FALSE, MZP=0.25, LC_F = "Auto", max_LC_F = 0.75, FG_Recov = 0.5, method="euclidean", VarPortion = 0.75,  ...) {
	#Data imputation function
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	if (MZP > 1 || MZP < 0) {
		stop("Invalid MZP.")
	}
	if (!is.logical(RowSamples)){
		stop("Invalid RowSamples.")
	}
	if (VarPortion >= 0.95 || VarPortion <= 0) {
		stop("Invalid VarPortion.")
	}
	if (LC_F > 0.95 || LC_F < 0.01) {
		if (LC_F != "Auto") {
			stop("Invalid LC_F.")
		}
	}
	datamatrix <- as.matrix(datamatrix)
	if (!RowSamples) {
		datamatrix <- t(datamatrix)
	}

	res.pca <- fast.prcomp(datamatrix)
	#Get number of pc with VarPortion
	num_pc2 <- 1
	Sumvar <- sum(res.pca$sdev^2)
	Sumwanted <- 0
	while(Sumwanted/Sumvar < VarPortion) {
		num_pc2 <- num_pc2 + 1
		Sumwanted <- sum(res.pca$sdev[1:num_pc2]^2, na.rm=TRUE)
	}
	DM <- as.matrix(Dist(res.pca$x[,1:num_pc2],method = method))
	DM[is.na(DM)] <- 1
	#prevent negatives
	if (min(DM) < 0) {
		DM <- DM + min(DM)
	}

	DM <- 1/(DM/rowSums(DM))
	DM[is.infinite(DM)] <- 1

	
	DM <- DM/rowSums(DM)

	#Filter dataset for genes that are suitable for DI
	datamatrix <- datamatrix[,colSums(datamatrix != 0) >= nrow(datamatrix) * MZP]
	if (LC_F == "Auto") {
		LC_F <- FindLCT_DI(datamatrix)
		if (max_LC_F  < LC_F) {
			LC_F <- max_LC_F
		}
		LC_F <- round(LC_F * (1 - FG_Recov),2)
	}
	MeanOrder <- order(colMeans(datamatrix),decreasing=FALSE)
	datamatrix <- datamatrix[,MeanOrder[floor(ncol(datamatrix) * LC_F + 1):ncol(datamatrix)]]
	

	Answer <- datamatrix
	for (i in 1:nrow(datamatrix)) {
		zeroes <- which(datamatrix[i,] == 0)
		yvec <- WcolMeans(datamatrix[,zeroes,drop=FALSE],DM[i,,drop=FALSE])
		Answer[i,zeroes] <- yvec
	}
	if (showinfo) {
		message("Number of PC chosen from PCA is ", num_pc2,".",appendLF=TRUE)
		message("Low expressing gene filtering threshold is set to ", LC_F, ".",appendLF=TRUE )
		flush.console()
	}
	if (!RowSamples) {
		Answer <- t(Answer)
	}
	return (Answer)
}
