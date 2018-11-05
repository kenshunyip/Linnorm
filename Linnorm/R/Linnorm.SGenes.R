#' Linnorm model stable gene selection tool.
#'
#' For datasets without spike-ins and for users who do not wish to rely on spike-ins, we provide this model stable gene selection tool.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets.
#' @param RowSamples	Logical. In the datamatrix, if each row is a sample and each row is a feature, set this to TRUE so that you don't need to transpose it. Linnorm works slightly faster with this argument set to TRUE, but it should be negligable for smaller datasets. Defaults to FALSE.
#' @param showinfo	Logical. Show algorithm running information. Defaults to FALSE.
#' @param minNonZeroPortion Double >=0, <= 1. Minimum non-Zero Portion Threshold. Genes not satisfying this threshold will be removed. For exmaple, if set to 0.3, genes without at least 30 percent of the samples being non-zero will be removed. Defaults to 0.3.
#' @param F_p	Double >=0, <= 1. Filter genes with standard deviation and skewness less than this p value. Defaults to 0.3173.
#' @param F_LC_Genes	Double >= 0.01, <= 0.95 or Character "Auto". Filter this portion of the lowest expressing genes. It can be determined automatically by setting to "Auto". Defaults to "Auto".
#' @param F_HC_Genes	Double >=0, <= 1. Filter this portion of the highest expressing genes. Defaults to 0.01.
#' @param max_F_LC	Double >=0, <= 0.95. When F_LC is set to auto, this is the maximum threshold that Linnorm would assign. Defaults to 0.75.
#' @details  This function selects stable genes from the dataset using the Linnorm's algorithm.
#' @return This function returns a data matrix that contains stable genes only.
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM Filter Stable
#' @export
#' @examples
#' #Obtain example matrix:
#' data(Islam2011)
#' #Transformation:
#' StableGenes <- Linnorm.SGenes(Islam2011)
#' @import
#' Rcpp
#' RcppArmadillo
Linnorm.SGenes <- function (datamatrix, RowSamples = FALSE, showinfo=FALSE, minNonZeroPortion = 0.3, F_p = 0.3173, F_LC_Genes = "Auto", F_HC_Genes = 0.01, max_F_LC = 0.75) {
	#Model stable gene selection
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	
	#data checking
	datamatrix <- as.matrix(datamatrix)
	RN <- 0
	CN <- 0
	if (RowSamples) {
		RN <- rownames(datamatrix)
		CN <- colnames(datamatrix)
		datamatrix <- XPM(datamatrix)
	} else {
		CN <- rownames(datamatrix)
		RN <- colnames(datamatrix)
		datamatrix <- tXPM(datamatrix) 
	}
	colnames(datamatrix) <- CN
	rownames(datamatrix) <- RN
	
	if (length(datamatrix[,1]) < 3) {
		stop("Number of samples is less than 3.")
	}
	if (length(datamatrix[1,]) < 500) {
		stop("Number of features is too small.")
	}
	if (minNonZeroPortion > 1 || minNonZeroPortion < 0) {
		stop("Invalid minNonZeroPortion.")
	}
	if (F_p > 1 || F_p < 0) {
		stop("Invalid F_p.")
	}
	if (F_LC_Genes > 0.75 || F_LC_Genes < 0.01) {
		if (F_LC_Genes != "Auto") {
			stop("Invalid F_LC_Genes.")
		}
	}
	if (F_HC_Genes > 0.75 || F_HC_Genes < 0.01) {
		stop("Invalid F_HC_Genes.")
	}
	if (anyNA(datamatrix)) {
		stop("Dataset contains NA.")
	}
	if (sum(which(datamatrix < 0)) != 0) {
		stop("Dataset contains negative number.")
	}
	if (max_F_LC > 0.95 || max_F_LC < 0) {
		stop("Invalid max_F_LC.")
	}
	if (!is.logical(RowSamples)){
		stop("Invalid RowSamples.")
	}
	if (!is.logical(showinfo)){
		stop("Invalid showinfo.")
	}
	#Find maxBound
	TheMean <- NZcolMeans(datamatrix)
	MeanOrder <- order(TheMean, decreasing = FALSE)
	numZero <- sum(TheMean == 0)
	fivepercent <- floor(0.05 * ncol(datamatrix)) + 1
	nonZero <- datamatrix[,MeanOrder[numZero:(numZero +fivepercent)]][which(datamatrix[,MeanOrder[numZero:(numZero +fivepercent)]] != 0)]
	maxBound <- length(nonZero)/sum(nonZero)
	
	#Get filter low count genes threhsold
	Keep <- 0
	if (nrow(datamatrix) * minNonZeroPortion < 3) {
		Keep <- which(colSums(datamatrix != 0) >= 3)
	} else {
		if (minNonZeroPortion == 0 || minNonZeroPortion == 1) {
			Keep <- which(colSums(datamatrix != 0) >= nrow(datamatrix) * minNonZeroPortion)
		} else {
			Keep <- which(colSums(datamatrix != 0) > nrow(datamatrix) * minNonZeroPortion)
		}
	}
	if (length(Keep) < 200) {
		stop("Given the current minNonZeroPortion threshold, the number of remaining feature (less than 200) is too small.")
	}
	LC_Threshold <- 0
	if (F_LC_Genes == "Auto") {
		LC_Threshold <- FindLCT(datamatrix[,Keep], maxBound)
		if (LC_Threshold > max_F_LC) {
			if (showinfo) {
				message(paste("Filter low count gene threshold is ", LC_Threshold, ". It is larger than max_F_LC, ", max_F_LC, ", which is now used.", sep=""))
			}
			LC_Threshold <- max_F_LC
		}
		F_LC_Genes <- LC_Threshold
		
		if (showinfo) {
			message(paste("Filter low count genes threshold is set to ", LC_Threshold, sep=""),appendLF=TRUE)
		}
	}
	if (F_LC_Genes + F_HC_Genes > 0.95){
		F_HC_Genes <- 0.01
		if (showinfo) {
			message(paste("F_HC_Genes Reset to ", F_HC_Genes, sep=""),appendLF=TRUE)
		}
	}
	
	#Filter dataset based on stdev and skewness
	datamatrix <- FirstFilter(datamatrix, minNonZeroPortion, L_F_p = F_p, L_F_LC_Genes = F_LC_Genes, L_F_HC_Genes = F_HC_Genes, spikein = NULL)
	
	if (!RowSamples) {
		datamatrix <- t(datamatrix)
	}
	return(datamatrix)
}

