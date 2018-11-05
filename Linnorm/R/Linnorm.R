#' Linnorm Normalizing Transformation Function
#'
#' This function performs the Linear model and normality based transformation method (Linnorm) for (sc)RNA-seq expression data or large scale count data.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets.
#' @param RowSamples	Logical. In the datamatrix, if each row is a sample and each row is a feature, set this to TRUE so that you don't need to transpose it. Linnorm works slightly faster with this argument set to TRUE, but it should be negligable for smaller datasets. Defaults to FALSE.
#' @param spikein	character vector. Names of the spike-in genes in the datamatrix. Defaults to NULL.
#' @param spikein_log2FC	Numeric vector. Log 2 fold change of the spike-in genes. Defaults to NULL.
#' @param showinfo	Logical. Show algorithm running information. Defaults to FALSE.
#' @param perturbation	Integer >=2 or "Auto". To search for an optimal minimal deviation parameter (please see the article), Linnorm uses the iterated local search algorithm which perturbs away from the initial local minimum. The range of the area searched in each perturbation is exponentially increased as the area get further away from the initial local minimum, which is determined by their index. This range is calculated by 10 * (perturbation ^ index). Defaults to "Auto".
#' @param Filter	Logical. Should Linnorm filter the dataset in the end results? Defaults to FALSE.
#' @param minNonZeroPortion Double >= 0.01, <= 0.95. Minimum non-Zero Portion Threshold. Genes not satisfying this threshold will be removed. For exmaple, if set to 0.3, genes without at least 30 percent of the samples being non-zero will be removed. Defaults to 0.3.
#' @param L_F_p	Double >= 0, <= 1. Filter genes with standard deviation and skewness less than this p value before applying Linnorm algorithm. Defaults to 0.3173.
#' @param L_F_LC_Genes	Double >= 0.01, <= 0.95 or Character "Auto". Filter this portion of the lowest expressing genes before applying Linnorm algorithm. It can be determined automatically by setting to "Auto". Defaults to "Auto".
#' @param L_F_HC_Genes	Double >= 0.01, <= 0.95. Filter this portion of the highest expressing genes before applying Linnorm algorithm. Defaults to 0.01.
#' @param BE_F_p	Double >=0, <= 1. Filter genes with standard deviation and skewness less than this p value before applying Linnorm's batch effect normalization algorithm. Defaults to 0.3173.
#' @param BE_F_LC_Genes	Double >= 0.01, <= 0.95 or Character "Auto". Filter this portion of the lowest expressing genes before applying Linnorm's batch effect normalization algorithm. It can be determined automatically by setting to "Auto". Defaults to "Auto".
#' @param BE_F_HC_Genes	Double >= 0.01, <= 0.95. Filter this portion of the highest expressing genes before applying Linnorm's batch effect normalization algorithm. Defaults to 0.01.
#' @param BE_strength	Double >0, <= 1. Before Linnorm transformation, how strongly should Linnorm normalize batch effects? Defaults to 0.5.
#' @param max_F_LC	Double >=0, <= 0.95. When L_F_LC or B_F_LC is set to auto, this is the maximum threshold that Linnorm would assign. Defaults to 0.75.
#' @param DataImputation	Logical. Perform data imputation on the dataset after transformation. Defaults to FALSE.
#' @param ... place holder for any new arguments.
#' @details  This function normalizes and transforms the input dataset using the Linnorm algorithm.
#' @return This function returns a transformed data matrix.
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric
#' @export
#' @examples
#' #Obtain example matrix:
#' data(LIHC)
#' #Transformation:
#' transformedExp <- Linnorm(LIHC)
#' @import
#' Rcpp
#' RcppArmadillo
Linnorm <- function(datamatrix, RowSamples = FALSE, spikein = NULL, spikein_log2FC = NULL, showinfo = FALSE, perturbation="Auto", Filter=FALSE, minNonZeroPortion = 0.3, L_F_p = 0.3173, L_F_LC_Genes = "Auto", L_F_HC_Genes = 0.01, BE_F_p = 0.3173, BE_F_LC_Genes = "Auto", BE_F_HC_Genes = 0.01, BE_strength = 0.5, max_F_LC=0.75, DataImputation = FALSE, ...) {
	#Linnorm transformation
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	
	#data checking
	datamatrix <- as.matrix(datamatrix)
	#Step 1: Relative Expression
	#Turn it into relative expression
	#Note that expdata does not have colnames and rownames now
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
	if (length(spikein) != length(spikein_log2FC)) {
		if (length(spikein_log2FC) == 0) {
			spikein_log2FC <- rep(0, length(spikein))
		} else {
			stop("spikein length must be the same as spikein_log2FC.")
		}
	} else {
		keep <- which(spikein_log2FC == 0)
		spikein <- spikein[keep]
	}
	if (perturbation < 2 && perturbation != "Auto") {
		stop("Invalid perturbation.")
	}
	if (minNonZeroPortion > 1 || minNonZeroPortion < 0) {
		stop("Invalid minNonZeroPortion.")
	}
	if (L_F_p >1 || L_F_p < 0) {
		stop("Invalid L_F_p.")
	}
	if (L_F_LC_Genes > 0.95 || L_F_LC_Genes < 0.01) {
		if (L_F_LC_Genes != "Auto") {
			stop("Invalid L_F_LC_Genes.")
		}
	}
	if (L_F_HC_Genes > 0.95 || L_F_HC_Genes < 0.01) {
		stop("Invalid L_F_HC_Genes.")
	}
	if (BE_F_p > 1 || BE_F_p < 0) {
		stop("Invalid BE_F_p.")
	}
	if (BE_F_LC_Genes > 0.95 || BE_F_LC_Genes < 0.01) {
		if (BE_F_LC_Genes != "Auto") {
			stop("Invalid BE_F_LC_Genes.")
		}
	}
	if (BE_F_HC_Genes > 0.95 || BE_F_HC_Genes < 0.01) {
		stop("Invalid BE_F_HC_Genes.")
	}
	if (BE_strength > 1 | BE_strength <= 0) {
		stop("Invalid BE_strength.")
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
	if (!is.logical(showinfo)){
		stop("Invalid showinfo.")
	}
	if (!is.logical(RowSamples)){
		stop("Invalid RowSamples.")
	}
	if (!is.logical(Filter)){
		stop("Invalid Filter.")
	}
	if (!is.logical(DataImputation)){
		stop("Invalid DataImputation.")
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
	if (BE_F_LC_Genes == "Auto" || L_F_LC_Genes == "Auto") {
		LC_Threshold <- FindLCT(datamatrix[,Keep], maxBound)
		if (LC_Threshold > max_F_LC) {
			if (showinfo) {
				message(paste("Filter low count gene threshold is ", LC_Threshold, ". It is larger than max_F_LC, ", max_F_LC, ", which is now used.", sep=""))
			}
			LC_Threshold <- max_F_LC
		}
		if (BE_F_LC_Genes == "Auto") {
			BE_F_LC_Genes <- LC_Threshold
		}
		if (L_F_LC_Genes == "Auto") {
			L_F_LC_Genes <- LC_Threshold
		}
		if (showinfo) {
			message(paste("Filter low count genes threshold is set to ", LC_Threshold, sep=""),appendLF=TRUE)
		}
	}
	if (L_F_LC_Genes + L_F_HC_Genes > 0.95){
		L_F_HC_Genes <- 0.01
		if (showinfo) {
			message(paste("L_F_HC_Genes Reset to ", L_F_HC_Genes, sep=""),appendLF=TRUE)
		}
	}
	if (BE_F_LC_Genes + BE_F_HC_Genes > 0.95){
		BE_F_HC_Genes <- 0.01
		if (showinfo) {
			message(paste("BE_F_HC_Genes Reset to ", BE_F_HC_Genes, sep=""),appendLF=TRUE)
		}
	}
	#Filter dataset and calculate lambda
	FilteredData <- FirstFilter(datamatrix, minNonZeroPortion, L_F_p = L_F_p, L_F_LC_Genes = L_F_LC_Genes, L_F_HC_Genes = L_F_HC_Genes, spikein = spikein)
	if (perturbation =="Auto") {
		perturbation <- round(maxBound^(1/3))+1
	}
	lambda <- LocateLambda(FilteredData, perturbation, maxBound)
	#Normalization
	if (BE_strength > 0) {
		datamatrix <- BatchEffectLinnorm1(datamatrix * lambda, minNonZeroPortion, BE_F_LC_Genes = BE_F_LC_Genes, BE_F_HC_Genes = BE_F_HC_Genes, BE_F_p = BE_F_p, BE_strength = BE_strength, spikein=spikein)
		datamatrix <- log1p(datamatrix)
	} else {
		datamatrix <- log1p(datamatrix * lambda)
	}
	x <- list(...)
	if (!is.null(x[['Internal']])) {
		Filter = TRUE
		LC_Threshold <- round(LC_Threshold * (1 - x$FG_Recov),2)
		minNonZeroPortion <- x$MZP
	}
	if (Filter || DataImputation) {
		if (Filter) {
			datamatrix <- datamatrix[,order(NZcolMeans(datamatrix),decreasing=FALSE)]
			datamatrix <- datamatrix[,colSums(datamatrix != 0) >= nrow(datamatrix) * minNonZeroPortion]
			Start <- floor(ncol(datamatrix) * LC_Threshold + 1)
			End <- ncol(datamatrix)
			Keep <- Start:End
			datamatrix <- datamatrix[,Keep]
		}
		if (DataImputation) {
			datamatrix <- Linnorm.DataImput(datamatrix, RowSamples=TRUE,showinfo=showinfo, ...)
		}
	}
	if (showinfo) {
		message("Lambda is ", lambda,".",appendLF=TRUE)
		message("perturbation is ", perturbation,".",appendLF=TRUE)
		flush.console()
	}
	if (RowSamples) {
		return (datamatrix)
	} else {
		return (t(datamatrix))
	}
	
}
