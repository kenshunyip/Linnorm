#' Linnorm-Hvar pipeline for highly variable gene discovery.
#'
#' This function first performs Linnorm transformation on the dataset. Then, it will perform highly variable gene discovery.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets.
#' @param RowSamples	Logical. In the datamatrix, if each row is a sample and each row is a feature, set this to TRUE so that you don't need to transpose it. Linnorm works slightly faster with this argument set to TRUE, but it should be negligable for smaller datasets. Defaults to FALSE.
#' @param spikein	character vector. Names of the spike-in genes in the datamatrix. Defaults to NULL.
#' @param spikein_log2FC	Numeric vector. Log 2 fold change of the spike-in genes. Defaults to NULL.
#' @param log.p	Logical. Output p/q values in log scale. Defaults to FALSE.
#' @param sig.value	Character. "p" or "q". Use p or q value for highlighting significant genes. Defaults to "p".
#' @param sig	Double >0, <= 1. Significant level of p or q value for plotting. Defaults to 0.05.
#' @param MZP Double >=0, <= 1. Minimum non-Zero Portion Threshold for this function. Genes not satisfying this threshold will be removed from HVG anlaysis. For exmaple, if set to 0.3, genes without at least 30 percent of the samples being non-zero will be removed. Defaults to 0.25.
#' @param FG_Recov	Double >=0, <= 1. In the low count gene filtering algorithm, recover this portion of genes that are filtered. Defaults to 0.5.
#' @param plot.title	Character. The plot's title. Defaults to "Mean vs SD plot".
#' @param ... arguments that will be passed into Linnorm's transformation function.
#' @details  This function discovers highly variable gene in the dataset using Linnorm transformation.
#' @return This function will output a list with the following objects:
##' \itemize{
##'  \item{Results:}{ A matrix with the results.}
##'  \item{plot:}{ Mean vs Standard Deviation Plot which highlights significant genes.}
##'  \item{Linnorm:}{ Linnorm transformed and filtered data matrix.}
##' }
#' @return The Results matrix has the following columns:
##' \itemize{
##'  \item{Transformed.Avg.Exp:}{ Average expression of non-zero Linnorm transformed data.}
##'  \item{Transformed.SD:}{ Standard deviation of non-zero Linnorm transformed data.}
##'  \item{Normalized.Log2.SD.Fold.Change:}{ Normalized log2 fold change of the gene's standard deviation.}
##'  \item{p.value:}{ p value of the statistical test.}
##'  \item{q.value:}{ q value/false discovery rate/adjusted p value of the statistical test.}
##' }
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric variance highly variable
#' @export
#' @examples
#' data(Islam2011)
#' results <- Linnorm.HVar(Islam2011)

Linnorm.HVar <- function(datamatrix, RowSamples = FALSE, spikein=NULL, spikein_log2FC=NULL, log.p=FALSE, sig.value="p", sig=0.05, MZP=0.25, FG_Recov=0.5, plot.title="Mean vs SD plot", ...) {
	#Highly variable gene analysis with Linnorm transformed dataset
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	datamatrix <- as.matrix(datamatrix)
	if (sig <= 0 || sig > 1) {
		stop("Invalid sig value.")
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
	if ( sig.value != "p" &&  sig.value != "q") {
		stop("Invalid sig.value.")
	}
	if (!is.logical(RowSamples)){
		stop("Invalid RowSamples.")
	}
	if (!RowSamples) {
		datamatrix <- t(datamatrix)
	}
	#Linnorm transformation
	datamatrix <- Linnorm(datamatrix, spikein=spikein, spikein_log2FC=spikein_log2FC, Internal=TRUE, MZP=MZP, FG_Recov=FG_Recov, RowSamples = TRUE, ...)
	#Filter genes based on number of non-zero values
	datamatrix <- datamatrix[,colSums(datamatrix != 0) >= 3]
	#Check available number of spike in genes.
	spikein <- spikein[which(spikein %in% colnames(datamatrix))]
	if (length(spikein) != 0 && length(spikein) < 10) {
		spikein = NULL
		warning("Too many Spikein are filtered. They will not be utilized.")
	}
	#Get mean and SD
	MeanSD <- NZcolMeanSD_acc(datamatrix)
	datamean <- MeanSD[1,]
	dataSD <- MeanSD[2,]
	#Loess Fit
	logitit <- loessFit(dataSD,datamean, weights=exp(datamean))
	#In case of negative fit, where negative stdev is impossible in our case, change them into the smallest stdev in the dataset
	logitit$fitted[which(logitit$fitted <= 0)] <- min(logitit$fitted[which(logitit$fitted > 0)])
	
	#Obtain Stdev ratios to adjust for technical noise.
	SDRatio <- as.numeric(log(dataSD/logitit$fitted,2))
	#normalize SDRatio
	LR <- LinearRegression(datamean,SDRatio)
	Residual <- (SDRatio - (LR$coefficients[[2]] * datamean + LR$coefficients[[1]]))
	LR2 <- LinearRegression(datamean,abs(Residual))
	SDRatio <- SDRatio * (LR2$coefficients[[2]] * datamean[1] + LR2$coefficients[[1]])/(LR2$coefficients[[2]] * datamean + LR2$coefficients[[1]])
	LR <- LinearRegression(datamean,SDRatio)
	Residual <- (SDRatio - (LR$coefficients[[2]] * datamean + LR$coefficients[[1]]))
	#Calculate p values
	pvalues <- 0
	spikes <- 0
	#if spike in list is provided, we use spike ins.
	if (length(spikein) < 10) {
		#Remove outlier
		SDRatio2 <- SDRatio[!SDRatio %in% boxplot.stats(SDRatio)$out]
		TheMean <- mean(SDRatio2)
		tdeno <- sqrt(sum((SDRatio2 - TheMean)^2)/(length(SDRatio2) - 2))
		pvalues <- pt((SDRatio - TheMean)/tdeno, df = length(SDRatio2) - 2, lower.tail = FALSE, log.p = TRUE)
	} else {
		spikes <- which(colnames(datamatrix) %in% spikein)
		SDRatio2 <- SDRatio[spikes]
		TheMean <- mean(SDRatio2)
		tdeno <- sqrt(sum((SDRatio2 - TheMean)^2)/(length(SDRatio2) - 2))
		pvalues <- pt((SDRatio - TheMean)/tdeno, df = length(SDRatio2) - 2, lower.tail = FALSE, log.p = TRUE)
		
	}
	#Organize results for output
	epvalues <- exp(pvalues)
	qvalues <- p.adjust(epvalues,"BH")
	dataSD <- dataSD
	results <- matrix(ncol=5, nrow=length(SDRatio))
	colnames(results) <- c("Transformed.Avg.Exp", "Transformed.SD", "Normalized.Log2.SD.Fold.Change", "p.value", "q.value")
	rownames(results) <- colnames(datamatrix)
	results[,1] <- datamean
	results[,2] <- dataSD
	results[,3] <- SDRatio
	if (log.p) {
		results[,4] <- pvalues
		results[,5] <- log(qvalues)
	} else {
		results[,4] <- epvalues
		results[,5] <- qvalues
	}
	#Plot Mean vs stdev
	groups <- rep("non-sig", length(SDRatio))
	if (sig.value == "p") {
		groups[which(epvalues <= sig)] <- "Significant"
	}
	if (sig.value == "q") {
		groups[which(qvalues <= sig)] <- "Significant"
	}
	if (spikes != 0) {
		groups[spikes] <- "Spike in"
		myColors <- c("blue","green","red")
		names(myColors) <- levels(groups)
	} else {
		myColors <- c("blue","red")
		names(myColors) <- levels(groups)
	}
	
	plotdata <- data.frame(mean=datamean,SD=dataSD,group=groups)
	render_plot <- ggplot_build(ggplot(plotdata, aes(x=mean, y=SD, color=group)) + geom_point(size = 1) + scale_x_continuous("Transformed Mean") + scale_y_continuous("Transformed Standard Deviation") + scale_colour_manual(name = "Sig",values = myColors) + ggtitle(plot.title) + theme(aspect.ratio=3/4))
	if (!RowSamples) {
		datamatrix <- t(datamatrix)
	}
	#Results for output
	listing <- list(results, render_plot, datamatrix)
	result <- setNames(listing, c("Results", "plot", "Linnorm"))
	return (result)
}


