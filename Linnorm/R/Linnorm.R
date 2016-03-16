#' Linnorm Normalization Function
#'
#' This function performs the Linear model and normality based normalization method (Linnorm) for RNA-seq expression data or large scale count data.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Undefined values such as NA are not supported.
#' @param showinfo Logical. Show lambda value calculated. Defaults to FALSE.
#' @param method	"default" or "lambda" The program will output the transformed matrix if the method is "default". If the method is "lambda", the program will output a lambda value.
#' @param minZeroPortion Double >=0, <= 1. Featuress without at least this portion of non-zero values will not be used in the calculation of normalizing parameter. Defaults to 2/3.
#' @param perturbation Integer >=2. To search for an optimal minimal deviation parameter (please see the article), Linnorm uses the iterated local search algorithm which perturbs away from the initial local minimum. The range of the area searched in each perturbation is exponentially increased as the area get further away from the initial local minimum, which is determined by their index. This range is calculated by 10 * (perturbation ^ index).
#' @details  If method is default, Linnorm outputs a transformed expression matrix. For users who wish to work with lambda instead, the output is a single lambda value. Please note that users with the lambda value can obtain a normalized Linnorm dataset by: log1p(lambda * datamatrix). There is no need to rerun the program if a lambda is already calculated.
#' @return This function returns either a transformed data matrix or a lambda value.
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric
#' @export
#' @examples
#' #Obtain example matrix:
#' library(seqc)
#' SampleA <- ILM_aceview_gene_BGI[,grepl("A_",colnames(ILM_aceview_gene_BGI))]
#' rownames(SampleA) <- ILM_aceview_gene_BGI[,2]
#' #Extract a portion of the matrix for an example
#' expMatrix <- SampleA[,1:3]
#' normalizedExp <- Linnorm(expMatrix)
#' normalizedExp <- Linnorm(expMatrix, method = "lambda")
#' @import
#' Rcpp
#' RcppArmadillo
Linnorm <- function(datamatrix, showinfo = FALSE, method="default",perturbation=10, minZeroPortion=2/3) {
	#data checking
	expdata <- as.matrix(datamatrix)
	if (length(expdata[1,]) < 2) {
		stop("Number of samples is less than 2.")
	}
	if (length(expdata[,1]) < 10) {
		stop("Number of features is too small.")
	}
	if (perturbation < 2) {
		stop("perturbation is too small.")
	}
	if (method != "default" && method != "lambda") {
		stop("Invalid algorithm value.")
	}
	if (minZeroPortion >1 || minZeroPortion < 0) {
		stop("Invalid minZeroPortion.")
	}
	#Step 1: Relative Expression
	#Turn it into relative expression
	for (i in seq_along(expdata[1,])) {
		expdata[,i] <- expdata[,i]/sum(expdata[,i])
	}
	#Save the original dataset before trimming.
	datamatrix <- expdata
	expdata <- expdata[order(rowMeans(expdata)),]	
	if (method == "default") {		
		#trim outliers
		expdata <- expdata[rowSums(expdata != 0) >= (length(expdata[1,]) * minZeroPortion),]
		
		X <- LocateLambda(expdata,perturbation)
		if (showinfo) {
			message("Lambda is ", X,".",appendLF=TRUE)
			flush.console()
		}
		return( log1p(datamatrix * X) )
	}
	if (method == "lambda") {
		#trim outliers
		expdata <- expdata[rowSums(expdata != 0) >= (length(expdata[1,]) * minZeroPortion),]
		
		X <- LocateLambda(expdata,perturbation)
		return( X )
	}
}


#' Linnorm-limma pipeline for Differentially Expression Analysis
#'
#' This function first performs Linnorm normalization on the dataset. Then, it will perform limma for DEG analysis. Finally, it will correct fold change outputs from limma results, that will be wrong otherwise. Please cite both Linnorm and limma when you use this function for publications.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Undefined values such as NA are not supported.
#' @param design	A design matrix required for limma. Please see limma's documentation or our vignettes for more detail.
#' @param output Character. "DEResults" or "Both". Set to "DEResults" to output a matrix that contains Differential Expression Analysis Results. Set to "Both" to output a list that contains both Differential Expression Analysis Results and the transformed data matrix.
#' @param noINF	Logical. Prevent generating INF in the fold change column by using Linnorm's lambda and adding one. If it is set to FALSE, INF will be generated if one of the conditions has zero expression. Defaults to TRUE. 
#' @param showinfo Logical. Show lambda value calculated. Defaults to FALSE.
#' @param minZeroPortion Double >=0, <= 1. Featuress without at least this portion of non-zero values will not be used in the calculation of normalizing parameter. Defaults to 2/3.
#' @param perturbation Integer >=2. To search for an optimal minimal deviation parameter (please see the article), Linnorm uses the iterated local search algorithm which perturbs away from the initial local minimum. The range of the area searched in each perturbation is exponentially increased as the area get further away from the initial local minimum, which is determined by their index. This range is calculated by 10 * (perturbation ^ index).
#' @param robust Logical. In the eBayes function of Limma, run with robust setting with TRUE or FALSE. Defaults to TRUE.
#' @details  This function performs both Linnorm and limma for users who are interested in differential expression analysis. Please note that if you directly use a Linnorm Nomralized dataset with limma, the output fold change and average expression with be wrong. (p values and adj.pvalues will be fine.) This is because the voom-limma pipeline assumes input to be in raw counts. This function is written to fix this problem.
#' @return If output is set to "DEResults", this function will output a matrix with Differntial Expression Analysis Results with the following columns:
##' \itemize{
##'  \item{logFC:}{ Log 2 Fold Change}
##'  \item{XPM:}{ Average Expression. If input is raw count or CPM, this column has the CPM unit. If input is RPKM, FPKM or TPM, this column has the TPM unit.}
##'  \item{t:}{ moderated t-statistic}
##'  \item{P.Value:}{ p value}
##'  \item{adj.P.Val:}{ Adjusted p value. This is also called False Discovery Rate or q value.}
##'  \item{B:}{ log odds that the feature is differential}
##' }
#' @return If output is set to Both, this function will output a list with the following objects:
##' \itemize{
##'  \item{DEResults:}{ Differntial Expression Analysis Results as described above.}
##'  \item{TMatrix:}{ A Linnorm Normalized Expression Matrix.}
##' }
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric limma
#' @export
#' @examples
#' #Obtain example matrix:
#' library(seqc)
#' SampleA <- ILM_aceview_gene_BGI[,grepl("A_",colnames(ILM_aceview_gene_BGI))]
#' rownames(SampleA) <- ILM_aceview_gene_BGI[,2]
#' SampleB <- ILM_aceview_gene_BGI[,grepl("B_",colnames(ILM_aceview_gene_BGI))]
#' rownames(SampleB) <- ILM_aceview_gene_BGI[,2]
#' #Extract a portion of the matrix for an example
#' expMatrix <- cbind(SampleA[,1:3], SampleB[,1:3])
#' designmatrix <- c(1,1,1,2,2,2)
#' designmatrix <- model.matrix(~ 0+factor(designmatrix))
#' colnames(designmatrix) <- c("group1", "group2")
#' rownames(designmatrix) <- colnames(expMatrix)
#' 
#' #Example 1
#' DEGResults <- Linnorm.limma(expMatrix, designmatrix)
#' #Example 2
#' DEGResults <- Linnorm.limma(expMatrix, designmatrix, output="Both")
Linnorm.limma <- function(datamatrix, design=NULL, output="DEResults", noINF=TRUE, showinfo = FALSE, perturbation=10, minZeroPortion=2/3, robust=TRUE) {
	expdata <- as.matrix(datamatrix)
	if (length(expdata[1,]) < 2) {
		stop("Number of samples is less than 2.")
	}
	if (length(expdata[,1]) < 10) {
		stop("Number of features is too small.")
	}
	if (perturbation < 2) {
		stop("perturbation is too small.")
	}
	if (minZeroPortion >1 || minZeroPortion < 0) {
		stop("Invalid minZeroPortion.")
	}
	if (is.null(design)) {
		stop("design is null.")
	}

	
	#Step 1: Relative Expression
	#Turn it into relative expression
	for (i in seq_along(expdata[1,])) {
		expdata[,i] <- expdata[,i]/sum(expdata[,i])
	}
	#Save the original dataset before trimming.
	datamatrix <- expdata
	expdata <- expdata[order(rowMeans(expdata)),]
	#trim outliers
	expdata <- expdata[rowSums(expdata != 0) >= (length(expdata[1,]) * minZeroPortion),]
	
	X <- LocateLambda(expdata,perturbation)
	if (showinfo) {
		message("Lambda is ", X,".",appendLF=TRUE)
		flush.console()
	}
	
	expdata <- log1p(datamatrix * X)
	
	#adjust design for further analysis
	CN <- c()
	for (i in seq_along(design[1,])) {
		CN <- c(CN, paste("group",i,sep=""))
	}
	colnames(design) <- CN
	fit2 <- lmFit(expdata,design)
	if (length(design[1,]) == 2) {
		contrast.matrix <- makeContrasts("group1-group2", levels=design)
	} else if (length(design[1,]) == 3) {
		contrast.matrix <- makeContrasts("group1-group2","group1-group3","group2-group3", levels=design)
	} else if (length(design[1,]) == 4) {
		contrast.matrix <- makeContrasts("group1-group2","group1-group3","group1-group4","group2-group3","group2-group4","group3-group4", levels=design)
	} else if (length(design[1,]) == 5) {
		contrast.matrix <- makeContrasts("group1-group2","group1-group3","group1-group4","group1-group5","group2-group3","group2-group4","group2-group5","group3-group4","group3-group5","group4-group5", levels=design)
	} else if (length(design[1,]) == 6) {
		contrast.matrix <- makeContrasts("group1-group2","group1-group3","group1-group4","group1-group5","group1-group6","group2-group3","group2-group4","group2-group5","group2-group6","group3-group4","group3-group5","group3-group6","group4-group5","group4-group6","group5-group6", levels=design)
	} else {
		stop("Error: The number of columns in design matrix is larger than 6.")
	}
	fit2 <- contrasts.fit(fit2, contrast.matrix)
	fit2 <- eBayes(fit2,robust=robust)
	limmaResults <- topTable(fit2, number=length(fit2$p.value), adjust.method="BH")
	datamatrix <- datamatrix[rownames(limmaResults),]
	if (length(design[1,]) == 2) {
		set1 <- as.numeric(which(design[,1] == 1))
		set2 <- as.numeric(which(design[,1] != 1))
		limmaResults[,2] <- unlist(apply(datamatrix,1,mean)) * 1000000
		if (noINF) {
			datamatrix <- datamatrix * X + 1
			limmaResults[,1] <- unlist(apply(datamatrix,1,function(x){return(log(mean(x[set1])/mean(x[set2] ),2) )}))
		} else {
			limmaResults[,1] <- unlist(apply(datamatrix,1,function(x){return(log(mean(x[set1])/mean(x[set2]),2) )}))
		}		
		colnames(limmaResults)[2] <- "XPM"
	} else {
		limmaResults[,2] <- unlist(apply(datamatrix,1,mean)) * 1000000
		colnames(limmaResults)[2] <- "XPM"
		limmaResults <- limmaResults[,-c(1)]
	}
	if (output=="DEResults") {
		return (limmaResults)
	}
	if (output=="Both") {
		listing <- list(limmaResults, expdata)
		results <- setNames(listing, c("DEResults", "TMatrix"))
		return (results)
	}
		
	
}


#' This function simulates a RNA-seq dataset based on a given distribution.
#' @param thisdata Matrix:	The matrix or data frame that contains your dataset. Each row is a gene and each column is a replicate. Undefined values such as NA are not supported. This program assumes that all columns are replicates of the same sample.
#' @param distribution Character: Defaults to "Poisson". This parameter controls the output distribution of the simulated RNA-seq dataset. It can be one of "Gamma" (Gamma distribution), "Poisson" (Poisson distribution), "LogNorm" (Log Normal distribution) or "NB" (Negative Binomial distribution).
#' @param NumRep Integer: The number of replicates. This is half of the number of output samples. Defaults to 3.
#' @param NumDiff Integer: The number of Differentially Changed Features. Defaults to 5000.
#' @param NumFea Integer: The number of Total Features. Defaults to 20000.
#' @param showinfo Logical: should we show data information on the console? Defaults to FALSE.
#' @param MaxLibSizelog2FC Double: The maximum library size difference from the mean that is allowed, in terms of log 2 fold change. Set to 0 to prevent program from generating library size differences. Defaults to 0.5.
#' @return This function returns a list that contains a matrix of count data in integer raw count and a vector that shows which genes are differentially expressed. In the matrix, each row is a gene and each column is a replicate. The first NumRep (see parameter) of the columns belong to sample 1, and the last NumRep (see parameter) of the columns belong to sample 2. There will be NumFea (see parameter) number of rows.
#' @keywords RNA-seq Raw Count Expression Simulation Gamma distribution Simulate Poisson "Log Normal" "Negative Binomial"
#' @export
#' @examples
#' #Obtain example matrix:
#' library(seqc)
#' SampleA <- ILM_aceview_gene_BGI[,grepl("A_",colnames(ILM_aceview_gene_BGI))]
#' rownames(SampleA) <- ILM_aceview_gene_BGI[,2]
#' #Extract a portion of the matrix for an example
#' expMatrix <- SampleA[,1:10]
#' simulateddata <- RnaXSim(expMatrix, distribution="NB", NumRep=3, NumDiff = 500, NumFea = 3000)
RnaXSim <- function(thisdata, distribution="Poisson", NumRep=3, NumDiff = 5000, NumFea = 20000, showinfo=FALSE, MaxLibSizelog2FC=0.5) {
	thisdata <- na.omit(as.matrix(thisdata))
	if (distribution != "Gamma" && distribution != "Poisson" && distribution != "LogNorm" && distribution != "NB") {
		stop("Invalid distribution.")
	}
	if (NumDiff < 0) {
		stop("Invalid NumDiff value.")
	}
	if (NumRep < 2) {
		stop("Invalid NumRep value.")
	}
	if (NumFea < 0) {
		stop("Invalid NumFea value.")
	}
	#Turn it into relative expression
	LibSize <- colSums(thisdata)
	for (i in seq_along(thisdata[1,])) {
		thisdata[,i] <- (thisdata[,i] * min(LibSize))/sum(thisdata[,i])
	}
	#sort and remove features with zeros and less than one.
	thisdata <- thisdata[rowMeans(thisdata) >= 1,]
	thisdata_ori <- thisdata
	
	if (distribution == "Gamma") {
		thisdata <- thisdata[order(rowMeans(thisdata)),]
		thisdata <- thisdata[rowSums(thisdata != 0) == length(thisdata[1,]),]
		meanList <- rowMeans(thisdata)
		b <- log(meanList)
		
		Klist <- vector(mode="numeric",length(thisdata[,1]))	
		for (i in seq_along(thisdata[,1])) {
			Klist[i] <- gammaShape(thisdata[i,])
		}
		
		if (length(Klist) == length(meanList)) {
			Fit <- lm((log(Klist))~(log(meanList)))
		} else {
			Klist <- vector(mode="numeric",length(thisdata[,1]))
			meanList <- vector(mode="numeric",length(thisdata[,1]))
			for (i in 1:length(thisdata[,1])) {
				Klist[i] <- gammaShape(thisdata[i,])
				meanList[i] <- mean(thisdata[i,])
			}
			Fit <- lm((log(Klist))~(log(meanList)))
		}
		constant <- Fit$coefficients[[1]]
		slope <- Fit$coefficients[[2]]
		Klist <- function (x) {
			return(exp(slope * log(x) + constant))
		}
		gammamatrix <- matrix(0, ncol=(2 * NumRep), nrow=NumFea)
		RN <- vector(mode="character", NumFea)
		for (i in 1:NumFea) {
			RN[i] <- paste("Gene",i,sep="")
		}
		rownames(gammamatrix) <- RN	

		Proportion <- NumDiff/NumFea
		if (Proportion > 0.9){
			stop("Error: NumDiff is too large. Proportion of Differential Feature is larger than 90%.")
		}
		if (Proportion <= 0){
			stop("Error: NumDiff is too small. Proportion of Differential Feature is smaller than 5%.")
		}
		#Minimum Fold change (FC) for pvalue to reach 0.05
		pvalue <- 1
		minBound <- 1
		#FC of 5 is sure to be significant, so there is no need to search for boundary. If it ever happens to be not so, it also serves as safety for the program, since FC > 5 not being significant doesn't make sense.
		maxBound <- 5
		midBound <- round((minBound + maxBound)/2,4)
		s <- median(rowMeans(thisdata_ori))
		theK <- Klist(s)
		theta1 <- s/theK
		NR2 <- NumRep
		design <- matrix(nrow=(NR2 * 2), ncol=2)
		colnames(design) <- c("SampleInfo", "Gam")
		col1 <- c()
		for (i in 1:NR2) {
			col1 <- c(col1, "Normal")
		}
		for (i in 1:NR2) {
			col1 <- c(col1, "Tumor")
		}
		design[,1] <- col1
		design <- as.data.frame(design)
		pvstore <- vector(mode="numeric",100)
		while (midBound != minBound || midBound != maxBound) {
			s2 <- s * midBound
			theK2 <- Klist(s2)
			theta2 <- s2/theK2
			#Do the test 100 times for an average p value.
			for (i in 1:100) {
				a <- as.numeric(rgamma(NR2,shape=theK,scale=theta1))
				b <- as.numeric(rgamma(NR2,shape=theK2,scale=theta2))
				#WilTest <- wilcox.test(a,b,alternative = c("two.sided"))
				#pvstore[i] <- as.numeric(WilTest[3])
				expressionData <- c(log1p(a),log1p(b))
				design[,2] <- as.numeric(expressionData)
				AnovaSum <- summary(aov(Gam~SampleInfo,data=design))
				pvstore[i] <- AnovaSum[[1]][,5][1]
				if (is.na(pvstore[i])) {
					pvstore[i] <- 1
				}
			}
			if (mean(pvstore, na.rm=TRUE) < 0.05) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		SigFC <- log(midBound)
		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find the stddev.
		minBound <- 0
		#FC of 100 is sure to be significant, so there is no need to search for boundary
		maxBound <- 100
		midBound <- (minBound + maxBound)/2
		Probability <- vector(mode="double",NumFea)
		while (midBound != minBound && midBound != maxBound) {
			normmodel <- rnorm(NumFea,0,midBound)
			for (i in seq_along(normmodel)) {
				Probability[i] <- 1 - (2 * (pnorm(abs(normmodel[i]),0,SigFC/2, lower.tail = TRUE) - 0.5 ) )
			}
			Probability[is.na(Probability)] <- 1
			Probability <- as.numeric(p.adjust(Probability,method ="BH"))
			answer1 <- Proportion - length(Probability[Probability < 0.05])/NumFea
			if (answer1 < 0) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		FCstddev <- midBound
		if (showinfo) {
			message("Significant Log 2 Fold Change for p value is ", log(exp(SigFC),2),".",appendLF=TRUE)
			message("Final Standard Deviation of log 2 fold change in the dataset is ", log(exp(FCstddev),2),".",appendLF=TRUE)
			flush.console()
		}
		
		
		#Non differential features will have mean of 0 and stddev of SigFC
		changes <- sort(rnorm(NumFea,0,FCstddev))
		begin <- round(NumFea/2 - NumDiff/2,0)
		ending <- round(NumFea/2 + NumDiff/2,0)
		unchanged <- changes[begin:ending]
		changing <- c(changes[1:begin],changes[ending:length(changes)])

		
		unchanged <- exp(unchanged)
		changing <- exp(changing)
		tobechanged <- sample(1:NumFea,NumDiff)
		for (i in 1:NumFea) {
			#sample from thisdata
			TDsample <- sample(1:length(thisdata_ori[,1]),1)
			thisrow <- thisdata_ori[TDsample,]
			dmean <- mean(thisrow)
			theK <- Klist(dmean)
			theta <- dmean/theK
			if (dmean == 0) {
				gammamatrix[i,] <- rep(0,length(gammamatrix[1,]))
			} else {
				if ( i %in% tobechanged) {
					if (sample(1:2,1) == 1) {
						gammamatrix[i,1:NumRep] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
						TBC <- sample(changing,1)
						dmean <- dmean * TBC
						theK <- Klist(dmean)
						theta <- dmean/theK
						gammamatrix[i,(NumRep + 1):length(gammamatrix[1,])] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
					} else {
						gammamatrix[i,(NumRep + 1):length(gammamatrix[1,])] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
						TBC <- sample(changing,1)
						dmean <- dmean * TBC
						theK <- Klist(dmean)
						theta <- dmean/theK
						gammamatrix[i,1:NumRep] <- sample(as.vector(rgamma(NumRep * 2, shape=theK,scale=theta)),NumRep)
					}
				} else {
					if (sample(1:2,1) == 1) {
						gammamatrix[i,1:NumRep] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
						TBC <- sample(unchanged,1)
						dmean <- dmean * TBC
						theK <- Klist(dmean)
						theta <- dmean/theK
						gammamatrix[i,(NumRep + 1):length(gammamatrix[1,])] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
					} else {
						gammamatrix[i,(NumRep + 1):length(gammamatrix[1,])] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
						TBC <- sample(unchanged,1)
						dmean <- dmean * TBC
						theK <- Klist(dmean)
						theta <- dmean/theK
						gammamatrix[i,1:NumRep] <- sample(as.vector(rgamma(NumRep * 2,shape=theK,scale=theta)),NumRep)
					}
				}
			}
		}
		gammamatrix[is.na(gammamatrix)] = 0
		gammamatrix[is.infinite(gammamatrix)] = 0
		gammamatrix[gammamatrix < 0] = 0
		changes <- 2^(seq(-MaxLibSizelog2FC,0,0.001) )
		thesechanges <- sample(changes,length(gammamatrix[1,]),replace=TRUE)

		for (i in seq_along(gammamatrix[1,])) {
			gammamatrix[,i] <- gammamatrix[,i] * thesechanges[i]
		}

		gammamatrix <-  as.matrix(floor(gammamatrix))
		listing <- list(gammamatrix,as.numeric(sort(tobechanged)))
		results <- setNames(listing, c("data", "DiffList"))
		return (results)
	}

	if (distribution == "Poisson") {			
		poismatrix <- matrix(0, ncol=(2 * NumRep), nrow=NumFea)
		RN <- vector(mode="character", NumFea)
		for (i in 1:NumFea) {
			RN[i] <- paste("Gene",i,sep="")
		}
		rownames(poismatrix) <- RN	

		Proportion <- NumDiff/NumFea
		if (Proportion > 0.9){
			stop("Error: NumDiff is too large. Proportion of Differential Feature is larger than 90%.")
		}
		if (Proportion <= 0){
			stop("Error: NumDiff is too small. Proportion of Differential Feature is smaller than 5%.")
		}
		#Minimum FC for pvalue to reach 0.05
		pvalue <- 1
		minBound <- 1
		#FC of 5 is sure to be significant, so there is no need to search for boundary. If it ever happens to be not so, it also serves as safety for the program, since FC > 5 not being significant doesn't make sense.
		maxBound <- 5
		midBound <- round((minBound + maxBound)/2,4)
		s <- median(rowMeans(thisdata))
		NR2 <- NumRep
		design <- matrix(nrow=(NR2 * 2), ncol=2)
		colnames(design) <- c("SampleInfo", "POIS")
		col1 <- c()
		for (i in 1:NR2) {
			col1 <- c(col1, "Normal")
		}
		for (i in 1:NR2) {
			col1 <- c(col1, "Tumor")
		}
		#####################################
		
		design[,1] <- col1
		design <- as.data.frame(design)
		pvstore <- vector(mode="numeric",100)
		while (midBound != minBound || midBound != maxBound) {
			s2 <- (s*midBound)
			#Do the test 100 times for an average p value.
			for (i in 1:100) {
				expressionData <- c(log1p(as.numeric(rpois(NR2,s))),log1p(as.numeric(rpois(NR2,s2))))
				design[,2] <- as.numeric(expressionData)
				AnovaSum <- summary(aov(POIS~SampleInfo,data=design))
				pvstore[i] <- AnovaSum[[1]][,5][1]
				#WilTest <- wilcox.test(as.numeric(rpois(NR2,s)),as.numeric(rpois(NR2,s2)),alternative = c("two.sided"))
				#pvstore[i] <- as.numeric(WilTest[3])
				if (is.na(pvstore[i])) {
					pvstore[i] <- 1
				}
			}
			if (mean(pvstore) < 0.05) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		SigFC <- log(midBound)

		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find the stddev.
		minBound <- 0
		#FC of 100 is sure to be significant, so there is no need to search for boundary
		maxBound <- 100
		midBound <- (minBound + maxBound)/2
		Probability <- vector(mode="double",NumFea)
		while (midBound != minBound && midBound != maxBound) {
			normmodel <- rnorm(NumFea,0,midBound)
			for (i in seq_along(normmodel)) {
				Probability[i] <- 1 - (2 * (pnorm(abs(normmodel[i]),0,SigFC/2, lower.tail = TRUE) - 0.5 ) )
			}
			Probability[is.na(Probability)] <- 1
			Probability <- as.numeric(p.adjust(Probability,method ="BH"))
			answer1 <- Proportion - length(Probability[Probability < 0.05])/NumFea
			if (answer1 < 0) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		FCstddev <- midBound
		if (showinfo) {
			message("Significant Log 2 Fold Change for p value is ", log(exp(SigFC),2),".",appendLF=TRUE)
			message("Final Standard Deviation of log 2 fold change in the dataset is ", log(exp(FCstddev),2),".",appendLF=TRUE)
			flush.console()
		}
		
		#Non differential features will have mean of 0 and stddev of SigFC
		changes <- sort(rnorm(NumFea,0,FCstddev))
		begin <- round(NumFea/2 - NumDiff/2,0)
		ending <- round(NumFea/2 + NumDiff/2,0)
		unchanged <- changes[begin:ending]
		changing <- c(changes[1:begin],changes[ending:length(changes)])

		
		unchanged <- exp(unchanged)
		changing <- exp(changing)
		tobechanged <- sample(1:NumFea,NumDiff)
		for (i in 1:NumFea) {
			#sample from thisdata
			TDsample <- sample(1:length(thisdata[,1]),1)
			dmean <- mean(as.numeric(thisdata[TDsample,]))
			if (dmean == 0) {
				poismatrix[i,] <- rep(0,length(poismatrix[1,]))
			} else {
				if ( i %in% tobechanged) {
					if (sample(1:2,1) == 1) {
						poismatrix[i,1:NumRep] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
						TBC <- sample(changing,1)
						dmean <- dmean * TBC
						poismatrix[i,(NumRep + 1):length(poismatrix[1,])] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
					} else {
						poismatrix[i,(NumRep + 1):length(poismatrix[1,])] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
						TBC <- sample(changing,1)
						dmean <- dmean * TBC
						poismatrix[i,1:NumRep] <- sample(as.vector(rpois(NumRep * 2, dmean)),NumRep)
					}
				} else {
					if (sample(1:2,1) == 1) {
						poismatrix[i,1:NumRep] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
						TBC <- sample(unchanged,1)
						dmean <- dmean * TBC
						poismatrix[i,(NumRep + 1):length(poismatrix[1,])] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
					} else {
						poismatrix[i,(NumRep + 1):length(poismatrix[1,])] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
						TBC <- sample(unchanged,1)
						dmean <- dmean * TBC
						poismatrix[i,1:NumRep] <- sample(as.vector(rpois(NumRep * 2, dmean)),NumRep)
					}
				}
			}
		}
		poismatrix[is.na(poismatrix)] = 0
		poismatrix[is.infinite(poismatrix)] = 0
		poismatrix[poismatrix < 0] = 0
		changes <- 2^(seq(-MaxLibSizelog2FC,0,0.001) )
		thesechanges <- sample(changes,length(poismatrix[1,]),replace=TRUE)
		
		for (i in seq_along(poismatrix[1,])) {
			poismatrix[,i] <- poismatrix[,i] * thesechanges[i]
		}

		poismatrix <-  as.matrix(floor(poismatrix))
		listing <- list(poismatrix,as.numeric(sort(tobechanged)))
		results <- setNames(listing, c("data", "DiffList"))

		return (results)
	}

	if (distribution == "LogNorm") {
		thisdata <- thisdata[order(rowMeans(thisdata)),]
		thisdata <- thisdata[rowSums(thisdata != 0) == length(thisdata[1,]),]
		
		
		Fit <- lm((log(rowSDs(thisdata)))~(log(rowMeans(thisdata))))
		MeanList <- rowMeans(thisdata)
		FindSD <- function(inputmean) {
			return (exp(log(inputmean) * Fit$coefficients[[2]] + Fit$coefficients[[1]] ))
		}
		lnormmatrix <- matrix(0, ncol=(2 * NumRep), nrow=NumFea)
		RN <- vector(mode="character", NumFea)
		for (i in 1:NumFea) {
			RN[i] <- paste("Gene",i,sep="")
		}
		rownames(lnormmatrix) <- RN	

		Proportion <- NumDiff/NumFea
		if (Proportion > 0.9){
			stop("Error: NumDiff is too large. Proportion of Differential Feature is larger than 90%.")
		}
		if (Proportion <= 0){
			stop("Error: NumDiff is too small. Proportion of Differential Feature is smaller than 5%.")
		}
		#Minimum FC for pvalue to reach 0.05
		pvalue <- 1
		minBound <- 1
		#FC of 5 is sure to be significant, so there is no need to search for boundary. If it ever happens to be not so, it also serves as safety for the program, since FC > 5 not being significant doesn't make sense.
		maxBound <- 5
		midBound <- round((minBound + maxBound)/2,4)
		mediandata <- thisdata[round(length(thisdata[,1])/2,0),]
		theMean <- median(rowMeans(thisdata_ori))
		theSD <- FindSD(theMean)
		LNMean <- log(theMean) - 0.5 * log((theSD/theMean)^2 + 1)
		LNSD <- (log((theSD/theMean)^2 + 1) )^0.5
		NR2 <- NumRep
		design <- matrix(nrow=(NR2 * 2), ncol=2)
		colnames(design) <- c("SampleInfo", "lnorm")
		col1 <- c()
		for (i in 1:NR2) {
			col1 <- c(col1, "Normal")
		}
		for (i in 1:NR2) {
			col1 <- c(col1, "Tumor")
		}
		design[,1] <- col1
		design <- as.data.frame(design)
		pvstore <- vector(mode="numeric",100)
		while (midBound != minBound || midBound != maxBound) {
			theMean2 <- theMean * midBound
			theSD2 <- FindSD(theMean2)
			LNMean2 <- log(theMean2) - 0.5 * log((theSD2/theMean2)^2 + 1)
			LNSD2 <- (log((theSD2/theMean2)^2 + 1) )^0.5
			#Do the test 100 times for an average p value.
			for (i in 1:100) {
				a <- abs(as.numeric(rlnorm(NR2,meanlog = LNMean, sdlog = LNSD)))
				b <- abs(as.numeric(rlnorm(NR2,meanlog = LNMean2, sdlog = LNSD2)))
				#WilTest <- wilcox.test(a,b,alternative = c("two.sided"))
				#pvstore[i] <- as.numeric(WilTest[3])
				expressionData <- c(log1p(a),log1p(b))
				design[,2] <- as.numeric(expressionData)
				AnovaSum <- summary(aov(lnorm~SampleInfo,data=design))
				pvstore[i] <- AnovaSum[[1]][,5][1]
				if (is.na(pvstore[i])) {
					pvstore[i] <- 1
				}
			}
			if (mean(pvstore) < 0.05) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		SigFC <- log(midBound)

		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find the stddev.
		minBound <- 0
		#FC of 100 is sure to be significant, so there is no need to search for boundary
		maxBound <- 100
		midBound <- (minBound + maxBound)/2
		Probability <- vector(mode="double",NumFea)
		while (midBound != minBound && midBound != maxBound) {
			normmodel <- rnorm(NumFea,0,midBound)
			for (i in seq_along(normmodel)) {
				Probability[i] <- 1 - (2 * (pnorm(abs(normmodel[i]),0,SigFC/2, lower.tail = TRUE) - 0.5 ) )
			}
			Probability[is.na(Probability)] <- 1
			Probability <- as.numeric(p.adjust(Probability,method ="BH"))
			answer1 <- Proportion - length(Probability[Probability < 0.05])/NumFea
			if (answer1 < 0) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		FCstddev <- midBound
		if (showinfo) {
			message("Significant Log 2 Fold Change for p value is ", log(exp(SigFC),2),".",appendLF=TRUE)
			message("Final Standard Deviation of log 2 fold change in the dataset is ", log(exp(FCstddev),2),".",appendLF=TRUE)
			flush.console()
		}
		
		#Non differential features will have mean of 0 and stddev of SigFC
		changes <- sort(rnorm(NumFea,0,FCstddev))
		begin <- round(NumFea/2 - NumDiff/2,0)
		ending <- round(NumFea/2 + NumDiff/2,0)
		unchanged <- changes[begin:ending]
		changing <- c(changes[1:begin],changes[ending:length(changes)])

		unchanged <- exp(unchanged)
		changing <- exp(changing)
		tobechanged <- sample(1:NumFea,NumDiff)
		for (i in 1:NumFea) {
			#sample from thisdata
			TDsample <- sample(1:length(thisdata_ori[,1]),1)
			thisrow <- thisdata_ori[TDsample,]
			theMean <- mean(thisrow)
			theSD <- FindSD(theMean)
			LNMean <- log(theMean) - 0.5 * log((theSD/theMean)^2 + 1)
			LNSD <- (log((theSD/theMean)^2 + 1) )^0.5
			if (theMean == 0) {
				lnormmatrix[i,] <- rep(0,length(lnormmatrix[1,]))
			} else {
				if ( i %in% tobechanged) {
					if (sample(1:2,1) == 1) {
						lnormmatrix[i,1:NumRep] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
						TBC <- sample(changing,1)
						theMean <- theMean * TBC
						theSD <- FindSD(theMean)
						LNMean <- log(theMean) - 0.5 * log((theSD/theMean)^2 + 1)
						LNSD <- (log((theSD/theMean)^2 + 1) )^0.5
						lnormmatrix[i,(NumRep + 1):length(lnormmatrix[1,])] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
					} else {
						lnormmatrix[i,(NumRep + 1):length(lnormmatrix[1,])] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
						TBC <- sample(changing,1)
						theMean <- theMean * TBC
						theSD <- FindSD(theMean)
						LNMean <- log(theMean) - 0.5 * log((theSD/theMean)^2 + 1)
						LNSD <- (log((theSD/theMean)^2 + 1) )^0.5
						lnormmatrix[i,1:NumRep] <- (sample(as.vector(rlnorm(NumRep * 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
					}
				} else {
					if (sample(1:2,1) == 1) {
						lnormmatrix[i,1:NumRep] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
						TBC <- sample(unchanged,1)
						theMean <- theMean * TBC
						theSD <- FindSD(theMean)
						LNMean <- log(theMean) - 0.5 * log((theSD/theMean)^2 + 1)
						LNSD <- (log((theSD/theMean)^2 + 1) )^0.5
						lnormmatrix[i,(NumRep + 1):length(lnormmatrix[1,])] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
					} else {
						lnormmatrix[i,(NumRep + 1):length(lnormmatrix[1,])] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
						TBC <- sample(unchanged,1)
						theMean <- theMean * TBC
						theSD <- FindSD(theMean)
						LNMean <- log(theMean) - 0.5 * log((theSD/theMean)^2 + 1)
						LNSD <- (log((theSD/theMean)^2 + 1) )^0.5
						lnormmatrix[i,1:NumRep] <- (sample(as.vector(rlnorm(NumRep * 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
					}
				}
			}
		}
		
		lnormmatrix[is.na(lnormmatrix)] = 0
		lnormmatrix[is.infinite(lnormmatrix)] = 0
		lnormmatrix[lnormmatrix < 0] = 0
		changes <- 2^(seq(-MaxLibSizelog2FC,0,0.001) )
		thesechanges <- sample(changes,length(lnormmatrix[1,]),replace=TRUE)
		for (i in seq_along(lnormmatrix[1,])) {
			lnormmatrix[,i] <- lnormmatrix[,i] * thesechanges[i]
		}
		lnormmatrix <-  as.matrix(floor(lnormmatrix))
		listing <- list(lnormmatrix,as.numeric(sort(tobechanged)))
		results <- setNames(listing, c("data", "DiffList"))
		return (results)
	}

	if (distribution == "NB") {
		thisdata_ori <- thisdata
		thisdata <- thisdata[order(rowMeans(thisdata)),]
		thisdata <- thisdata[rowSums(thisdata != 0) == length(thisdata[1,]),]
		
		MeanList <- vector(mode="numeric",length(thisdata[,1]))
		d <- vector(mode="numeric",length(thisdata[,1]))
		
		for (i in seq_along(thisdata[,1])){
			x <- thisdata[i,]
			MeanList[i] <- mean(x)
			d[i] <- as.numeric(unlist((glm.nb(x~1))[[24]]))
		}
		MeanList <- MeanList[!is.na(d)]
		d <- d[!is.na(d)]

		Fit <- lm((log(d))~(log(MeanList)))
		FindDispersion <- function(inputmean) {
			return (exp(log(inputmean) * Fit$coefficients[[2]] + Fit$coefficients[[1]] ))
		}
		
		nbmatrix <- matrix(0, ncol=(2 * NumRep), nrow=NumFea)
		RN <- vector(mode="character", NumFea)
		for (i in 1:NumFea) {
			RN[i] <- paste("Gene",i,sep="")
		}
		rownames(nbmatrix) <- RN	
		
		Proportion <- NumDiff/NumFea
		if (Proportion > 0.9){
			stop("Error: NumDiff is too large. Proportion of Differential Feature is larger than 90%.")
		}
		if (Proportion <= 0){
			stop("Error: NumDiff is too small. Proportion of Differential Feature is smaller than 5%.")
		}
		#Minimum FC for pvalue to reach 0.05
		pvalue <- 1
		minBound <- 1
		#FC of 5 is sure to be significant, so there is no need to search for boundary. If it ever happens to be not so, it also serves as safety for the program, since FC > 5 not being significant doesn't make sense.
		maxBound <- 5
		midBound <- round((minBound + maxBound)/2,4)
		themean <- median(rowMeans(thisdata_ori))
		theDis <- FindDispersion(themean)
		
		NR2 <- NumRep
		design <- matrix(nrow=(NR2 * 2), ncol=2)
		colnames(design) <- c("SampleInfo", "Gam")
		col1 <- c()
		for (i in 1:NR2) {
			col1 <- c(col1, "Normal")
		}
		for (i in 1:NR2) {
			col1 <- c(col1, "Tumor")
		}
		design[,1] <- col1
		design <- as.data.frame(design)
		pvstore <- vector(mode="numeric",100)
		while (midBound != minBound || midBound != maxBound) {
			themean2 <- themean * midBound
			theDis2 <- FindDispersion(themean2)
			
			#Do the test 100 times for an average p value.
			for (i in 1:100) {
				a <- as.numeric(rnbinom(NR2,size=theDis,mu=themean))
				b <- as.numeric(rnbinom(NR2,size=theDis2,mu=themean2))
				#WilTest <- wilcox.test(a,b,alternative = c("two.sided"))
				#pvstore[i] <- as.numeric(WilTest[3])
				expressionData <- c(log1p(a),log1p(b))
				design[,2] <- as.numeric(expressionData)
				AnovaSum <- summary(aov(Gam~SampleInfo,data=design))
				pvstore[i] <- AnovaSum[[1]][,5][1]
				if (is.na(pvstore[i])) {
					pvstore[i] <- 1
				}
			}
			if (mean(pvstore, na.rm=TRUE) < 0.05) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		SigFC <- log(midBound)

		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find the stddev.
		minBound <- 0
		#FC of 100 is sure to be significant, so there is no need to search for boundary
		maxBound <- 100
		midBound <- (minBound + maxBound)/2
		Probability <- vector(mode="double",NumFea)
		while (midBound != minBound && midBound != maxBound) {
			normmodel <- rnorm(NumFea,0,midBound)
			for (i in seq_along(normmodel)) {
				Probability[i] <- 1 - (2 * (pnorm(abs(normmodel[i]),0,SigFC/2, lower.tail = TRUE) - 0.5 ) )
			}
			Probability[is.na(Probability)] <- 1
			Probability <- as.numeric(p.adjust(Probability,method ="BH"))
			answer1 <- Proportion - length(Probability[Probability < 0.05])/NumFea
			if (answer1 < 0) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		FCstddev <- midBound
		if (showinfo) {
			message("Significant Log 2 Fold Change for p value is ", log(exp(SigFC),2),".",appendLF=TRUE)
			message("Final Standard Deviation of log 2 fold change in the dataset is ", log(exp(FCstddev),2),".",appendLF=TRUE)
			flush.console()
		}
		
		#Non differential features will have mean of 0 and stddev of SigFC
		changes <- sort(rnorm(NumFea,0,FCstddev))
		begin <- round(NumFea/2 - NumDiff/2,0)
		ending <- round(NumFea/2 + NumDiff/2,0)
		unchanged <- changes[begin:ending]
		changing <- c(changes[1:begin],changes[ending:length(changes)])

		
		unchanged <- exp(unchanged)
		changing <- exp(changing)
		tobechanged <- sample(1:NumFea,NumDiff)
		for (i in 1:NumFea) {
			#sample from thisdata
			TDsample <- sample(1:length(thisdata_ori[,1]),1)
			thisrow <- thisdata_ori[TDsample,]
			themean <- mean(thisrow)
			theDis <- FindDispersion(themean)
		
			
			if (themean == 0) {
				nbmatrix[i,] <- rep(0,length(nbmatrix[1,]))
			} else {
				if ( i %in% tobechanged) {
					if (sample(1:2,1) == 1) {
						nbmatrix[i,1:NumRep] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
						TBC <- sample(changing,1)
						themean <- themean * TBC
						theDis <- FindDispersion(themean)
						nbmatrix[i,(NumRep + 1):length(nbmatrix[1,])] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
					} else {
						nbmatrix[i,(NumRep + 1):length(nbmatrix[1,])] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
						TBC <- sample(changing,1)
						themean <- themean * TBC
						theDis <- FindDispersion(themean)
						nbmatrix[i,1:NumRep] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
					}
				} else {
					if (sample(1:2,1) == 1) {
						nbmatrix[i,1:NumRep] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
						TBC <- sample(unchanged,1)
						themean <- themean * TBC
						theDis <- FindDispersion(themean)
						nbmatrix[i,(NumRep + 1):length(nbmatrix[1,])] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
					} else {
						nbmatrix[i,(NumRep + 1):length(nbmatrix[1,])] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
						TBC <- sample(unchanged,1)
						themean <- themean * TBC
						theDis <- FindDispersion(themean)
						nbmatrix[i,1:NumRep] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
					}
				}
			}
		}
		nbmatrix[is.na(nbmatrix)] = 0
		nbmatrix[is.infinite(nbmatrix)] = 0
		nbmatrix[nbmatrix < 0] = 0
		changes <- 2^(seq(-MaxLibSizelog2FC,0,0.001) )
		thesechanges <- sample(changes,length(nbmatrix[1,]),replace=TRUE)

		for (i in seq_along(nbmatrix[1,])) {
			nbmatrix[,i] <- nbmatrix[,i] * thesechanges[i]
		}

		nbmatrix <-  as.matrix(floor(nbmatrix))
		listing <- list(nbmatrix,as.numeric(sort(tobechanged)))
		results <- setNames(listing, c("data", "DiffList"))

		return (results)
	}
	
}

