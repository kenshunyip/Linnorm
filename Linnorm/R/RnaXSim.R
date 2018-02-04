#' This function simulates an RNA-seq dataset based on a given distribution.
#' @param datamatrix	Matrix. The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets.This program assumes that all columns are replicates of the same sample.
#' @param distribution	Character: Defaults to "Poisson". This parameter controls the output distribution of the simulated RNA-seq dataset. It can be one of "Gamma" (Gamma distribution), "Poisson" (Poisson distribution), "LogNorm" (Log Normal distribution) or "NB" (Negative Binomial distribution).
#' @param NumRep	Integer: The number of replicates. This is half of the number of output samples. Defaults to 3.
#' @param NumDiff	Integer: The number of Differentially Changed Features. Defaults to 2000.
#' @param NumFea	Integer: The number of Total Features. Defaults to 20000.
#' @param showinfo	Logical: should we show data information on the console? Defaults to FALSE.
#' @param DEGlog2FC	"Auto" or Double: log 2 fold change threshold that defines differentially expressed genes. If set to "Auto," DEGlog2FC is defined at the level where ANOVA can get a q value of 0.05 with the average expression, where the data values are log1p transformed. Defaults to "Auto".
#' @param MaxLibSizelog2FC	Double: The maximum library size difference from the mean that is allowed, in terms of log 2 fold change. Set to 0 to prevent program from generating library size differences. Defaults to 0.5.
#' @return This function returns a list that contains a matrix of count data in integer raw count and a vector that shows which genes are differentially expressed. In the matrix, each row is a gene and each column is a replicate. The first NumRep (see parameter) of the columns belong to sample 1, and the last NumRep (see parameter) of the columns belong to sample 2. There will be NumFea (see parameter) number of rows. The top NumCorr of genes will be positively or negatively correlated with each other (randomly); and they are evenly separated into groups. Each group is not intended to be correlated to each other, but, by chance, it can happen.
#' @keywords RNA-seq Raw Count Expression Simulation Gamma distribution Simulate Poisson "Log Normal" "Negative Binomial"
#' @export
#' @examples
#' #Obtain example matrix:
#' data(SEQC)
#' expMatrix <- SEQC
#' #Example for Negative Binomial distribution
#' simulateddata <- RnaXSim(expMatrix, distribution="NB", NumRep=5, NumDiff = 200, NumFea = 2000)
#' #Example for Poisson distribution
#' simulateddata <- RnaXSim(expMatrix, distribution="Poisson", NumRep=5, NumDiff = 200, NumFea = 2000)
#' #Example for Log Normal distribution
#' simulateddata <- RnaXSim(expMatrix, distribution="LogNorm", NumRep=5, NumDiff = 200, NumFea = 2000)
#' #Example for Gamma distribution
#' simulateddata <- RnaXSim(expMatrix, distribution="Gamma", NumRep=5, NumDiff = 200, NumFea = 2000)
RnaXSim <- function(datamatrix, distribution="NB", NumRep=5, NumDiff = 2000, NumFea = 20000, showinfo=FALSE, DEGlog2FC="Auto", MaxLibSizelog2FC=0.5) {
	#RNA-seq data simulation for DEG analysis evaluation
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	datamatrix <- na.omit(as.matrix(datamatrix))
	if (length(datamatrix[1,]) < 3) {
		stop("Number of samples is less than 3.")
	}
	if (length(datamatrix[,1]) < 500) {
		stop("Number of features is too small.")
	}
	if (distribution != "Gamma" && distribution != "Poisson" && distribution != "LogNorm" && distribution != "NB") {
		stop("Invalid distribution.")
	}
	if (NumDiff < 0) {
		stop("Invalid NumDiff value.")
	}
	if (NumRep < 5) {
		stop("Invalid NumRep value.")
	}
	if (NumFea < 0) {
		stop("Invalid NumFea value.")
	}
	if (anyNA(datamatrix)) {
		stop("Dataset contains NA.")
	}
	if (sum(which(datamatrix < 0)) != 0) {
		stop("Dataset contains negative number.")
	}
	
	#Turn it into relative expression
	LibSize <- colSums(datamatrix)
	for (i in seq_along(datamatrix[1,])) {
		datamatrix[,i] <- (datamatrix[,i] * max(LibSize))/sum(datamatrix[,i])
	}
	#sort and remove features with zeros and less than one.
	#datamatrix <- datamatrix[rowMeans(datamatrix) >= 1,]
	datamatrix <- datamatrix[rowSums(datamatrix != 0) == length(datamatrix[1,]),]

	datamatrix <- datamatrix[order(rowMeans(datamatrix)),]
	
	if (distribution == "Gamma") {
		return (GammaSim(datamatrix, NumRep=NumRep, NumDiff = NumDiff, NumFea = NumFea, showinfo=showinfo, MaxLibSizelog2FC=MaxLibSizelog2FC, DEGlog2FC=DEGlog2FC))
	}

	if (distribution == "Poisson") {
		return (PoissonSim(datamatrix, NumRep=NumRep, NumDiff = NumDiff, NumFea = NumFea, showinfo=showinfo, MaxLibSizelog2FC=MaxLibSizelog2FC, DEGlog2FC=DEGlog2FC))
	}

	if (distribution == "LogNorm") {
		return (LogNormSim(datamatrix, NumRep=NumRep, NumDiff = NumDiff, NumFea = NumFea, showinfo=showinfo, MaxLibSizelog2FC=MaxLibSizelog2FC, DEGlog2FC=DEGlog2FC))
	}

	if (distribution == "NB") {
		return (NBSim(datamatrix, NumRep=NumRep, NumDiff = NumDiff, NumFea = NumFea, showinfo=showinfo, MaxLibSizelog2FC=MaxLibSizelog2FC, DEGlog2FC=DEGlog2FC))
	}
	
}

GammaSim <- function(thisdata_ori, NumRep=5, NumDiff = 2000, NumFea = 20000, showinfo=FALSE, MaxLibSizelog2FC=0.5, DEGlog2FC="Auto") {
	#RNA-seq data simulation for DEG analysis evaluation
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	#Capture distribution from the dataset for simulation
	#K is the shape parameter from Gamma distribution
	MeanSD <- rowMeanSD(thisdata_ori)
	Fit <- LinearRegression(log(MeanSD[1,]),log(MeanSD[2,]^2))
	FindVar <- function(inputmean) {
		return (exp(log(inputmean) * Fit$coefficients[[2]] + Fit$coefficients[[1]] ))
	}
	
	gammamatrix <- matrix(0, ncol=(2 * NumRep), nrow=NumFea)
	RN <- vector(mode="character", NumFea)
	for (i in 1:NumFea) {
		RN[i] <- paste("Gene",i,sep="")
	}
	rownames(gammamatrix) <- RN	
	
	
	#Obtain Significant fold change threshold
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
	#FC of 5 is sure to be significant, so there is no need to search for boundary. If it ever happens to be not so, it also serves as safety for the program, since FC > 5  not being significant doesn't make sense.
	maxBound <- 5
	midBound <- round((minBound + maxBound)/2,4)
	s <- median(MeanSD[1,])
	theta1 <- FindVar(s)/s
	theK <- s/theta1
	
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
	if (DEGlog2FC != "Auto") {
		SigFC <- log(2^DEGlog2FC)
		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find its stddev.
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
	} else { 
		while (midBound != minBound || midBound != maxBound) {
			s2 <- s * midBound
			theta2 <- FindVar(s2)/s2
			theK2 <- s2/theta2
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
		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find its stddev.
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
	}
	#Define Fold Change of Genes.
	#Non differential features will have mean of 0 and stddev of SigFC
	changes <- sort(rnorm(NumFea,0,FCstddev))
	begin <- round(NumDiff/2,0)
	ending <- round(NumFea - NumDiff/2,0)
	unchanged <- changes[begin:ending]
	changing <- c(changes[1:begin],changes[ending:length(changes)])

	Group1 <- c()
	
	unchanged <- exp(unchanged)
	changing <- exp(changing)
	tobechanged <- sample(1:NumFea,NumDiff)
	for (i in 1:NumFea) {
		#sample from thisdata
		TDsample <- sample(1:length(thisdata_ori[,1]),1)
		thisrow <- thisdata_ori[TDsample,]
		dmean <- mean(thisrow)
		theta <- FindVar(dmean)/dmean
		theK <- dmean/theta
		if (dmean == 0) {
			gammamatrix[i,] <- rep(0,length(gammamatrix[1,]))
		} else {
			if ( i %in% tobechanged) {
				if (sample(1:2,1) == 1) {
					gammamatrix[i,1:NumRep] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
					TBC <- sample(changing,1)
					if (TBC > 1) {
						Group1 <- c(Group1, i)
					}
					dmean <- dmean * TBC
					theta <- FindVar(dmean)/dmean
					theK <- dmean/theta
					gammamatrix[i,(NumRep + 1):length(gammamatrix[1,])] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
				} else {
					gammamatrix[i,(NumRep + 1):length(gammamatrix[1,])] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
					TBC <- sample(changing,1)
					if (TBC < 1) {
						Group1 <- c(Group1, i)
					}
					dmean <- dmean * TBC
					theta <- FindVar(dmean)/dmean
					theK <- dmean/theta
					gammamatrix[i,1:NumRep] <- sample(as.vector(rgamma(NumRep * 2, shape=theK,scale=theta)),NumRep)
				}
			} else {
				if (sample(1:2,1) == 1) {
					gammamatrix[i,1:NumRep] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
					TBC <- sample(unchanged,1)
					dmean <- dmean * TBC
					theta <- FindVar(dmean)/dmean
					theK <- dmean/theta
					gammamatrix[i,(NumRep + 1):length(gammamatrix[1,])] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
				} else {
					gammamatrix[i,(NumRep + 1):length(gammamatrix[1,])] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
					TBC <- sample(unchanged,1)
					dmean <- dmean * TBC
					theta <- FindVar(dmean)/dmean
					theK <- dmean/theta
					gammamatrix[i,1:NumRep] <- sample(as.vector(rgamma(NumRep * 2,shape=theK,scale=theta)),NumRep)
				}
			}
		}
	}
	gammamatrix[is.na(gammamatrix)] = 0
	gammamatrix[is.infinite(gammamatrix)] = 0
	gammamatrix[gammamatrix < 0] = 0
	
	
	#Induce library size differences
	changes <- 2^(seq(-MaxLibSizelog2FC,0,0.001) )
	thesechanges <- sample(changes,length(gammamatrix[1,]),replace=TRUE)

	for (i in seq_along(gammamatrix[1,])) {
		gammamatrix[,i] <- gammamatrix[,i] * thesechanges[i]
	}

	gammamatrix <-  as.matrix(floor(gammamatrix))
	colnames(gammamatrix) <- paste("Sample",1:ncol(gammamatrix),sep="_")
	listing <- list(gammamatrix,as.numeric(sort(tobechanged)), as.numeric(sort(Group1)), as.numeric(sort( tobechanged[!(tobechanged %in% Group1)])))
	results <- setNames(listing, c("data", "DiffList", "CorGroup1", "CorGroup2"))
	return (results)
}

PoissonSim <- function(thisdata, NumRep=5, NumDiff = 2000, NumFea = 20000, showinfo=FALSE, MaxLibSizelog2FC=0.5, DEGlog2FC="Auto") {
	#RNA-seq data simulation for DEG analysis evaluation
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	thisdata_ori <- thisdata
	poismatrix <- matrix(0, ncol=(2 * NumRep), nrow=NumFea)
	RN <- vector(mode="character", NumFea)
	for (i in 1:NumFea) {
		RN[i] <- paste("Gene",i,sep="")
	}
	rownames(poismatrix) <- RN	
	
	#Obtain Significant fold change threshold
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
	#FC of 5 is sure to be significant, so there is no need to search for boundary. If it ever happens to be not so, it also serves as safety for the program, since FC > 5  not being significant doesn't make sense.
	maxBound <- 5
	midBound <- round((minBound + maxBound)/2,4)
	s <- median(rowMeans(thisdata_ori))
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
	
	if (DEGlog2FC != "Auto") {
		SigFC <- log(2^DEGlog2FC)
		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find its stddev.
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
	} else { 
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
	}
	#Non differential features will have mean of 0 and stddev of SigFC
	changes <- sort(rnorm(NumFea,0,FCstddev))
	begin <- round(NumDiff/2,0)
	ending <- round(NumFea - NumDiff/2,0)
	unchanged <- changes[begin:ending]
	changing <- c(changes[1:begin],changes[ending:length(changes)])

	Group1 <- c()
	unchanged <- exp(unchanged)
	changing <- exp(changing)
	tobechanged <- sample(1:NumFea,NumDiff)
	for (i in 1:NumFea) {
		#sample from thisdata
		TDsample <- sample(1:length(thisdata_ori[,1]),1)
		dmean <- mean(as.numeric(thisdata_ori[TDsample,]))
		if (dmean == 0) {
			poismatrix[i,] <- rep(0,length(poismatrix[1,]))
		} else {
			if ( i %in% tobechanged) {
				if (sample(1:2,1) == 1) {
					poismatrix[i,1:NumRep] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
					TBC <- sample(changing,1)
					if (TBC > 1) {
						Group1 <- c(Group1, i)
					}
					dmean <- dmean * TBC
					poismatrix[i,(NumRep + 1):length(poismatrix[1,])] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
				} else {
					poismatrix[i,(NumRep + 1):length(poismatrix[1,])] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
					TBC <- sample(changing,1)
					if (TBC < 1) {
						Group1 <- c(Group1, i)
					}
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
	

	#Induce library size differences
	changes <- 2^(seq(-MaxLibSizelog2FC,0,0.001) )
	thesechanges <- sample(changes,length(poismatrix[1,]),replace=TRUE)
	
	for (i in seq_along(poismatrix[1,])) {
		poismatrix[,i] <- poismatrix[,i] * thesechanges[i]
	}

	poismatrix <-  as.matrix(floor(poismatrix))
	colnames(poismatrix) <- paste("Sample",1:ncol(poismatrix),sep="_")
	listing <- list(poismatrix,as.numeric(sort(tobechanged)), as.numeric(sort(Group1)), as.numeric(sort( tobechanged[!(tobechanged %in% Group1)])))
	results <- setNames(listing, c("data", "DiffList", "CorGroup1", "CorGroup2"))

	return (results)
}

LogNormSim <- function(thisdata_ori, NumRep=5, NumDiff = 2000, NumFea = 20000, showinfo=FALSE, MaxLibSizelog2FC=0.5, DEGlog2FC="Auto") {
	#RNA-seq data simulation for DEG analysis evaluation
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	#Capture distribution from the dataset for simulation
	MeanSD <- rowMeanSD(thisdata_ori)
	Fit <- LinearRegression(log(MeanSD[1,]),log(MeanSD[2,]))
	FindSD <- function(inputmean) {
		return (exp(log(inputmean) * Fit$coefficients[[2]] + Fit$coefficients[[1]] ))
	}
	lnormmatrix <- matrix(0, ncol=(2 * NumRep), nrow=NumFea)
	RN <- vector(mode="character", NumFea)
	for (i in 1:NumFea) {
		RN[i] <- paste("Gene",i,sep="")
	}
	rownames(lnormmatrix) <- RN	
	
	
	#Obtain Significant fold change threshold
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
	#FC of 5 is sure to be significant, so there is no need to search for boundary. If it ever happens to be not so, it also serves as safety for the program, since FC > 5  not being significant doesn't make sense.
	maxBound <- 5
	midBound <- round((minBound + maxBound)/2,4)
	theMean <- median(MeanSD[1,])
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
	if (DEGlog2FC != "Auto") {
		SigFC <- log(2^DEGlog2FC)
		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find its stddev.
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
	} else { 
		while (midBound != minBound || midBound != maxBound) {
			theMean2 <- theMean * midBound
			theSD2 <- FindSD(theMean2)
			LNMean2 <- log(theMean2) - 0.5 * log((theSD2/theMean2)^2 + 1)
			LNSD2 <- (log((theSD2/theMean2)^2 + 1) )^0.5
			#Do the test 100 times for an average p value.
			for (i in 1:100) {
				a <- as.numeric(rlnorm(NR2,meanlog = LNMean, sdlog = LNSD))
				b <- as.numeric(rlnorm(NR2,meanlog = LNMean2, sdlog = LNSD2))
				a[a < 0] <- 0
				b[b < 0] <- 0
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
	}
	#Non differential features will have mean of 0 and stddev of SigFC
	changes <- sort(rnorm(NumFea,0,FCstddev))
	begin <- round(NumDiff/2,0)
	ending <- round(NumFea - NumDiff/2,0)
	unchanged <- changes[begin:ending]
	changing <- c(changes[1:begin],changes[ending:length(changes)])
	Group1 <- c()
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
					if (TBC > 1) {
						Group1 <- c(Group1, i)
					}
					theMean <- theMean * TBC
					theSD <- FindSD(theMean)
					LNMean <- log(theMean) - 0.5 * log((theSD/theMean)^2 + 1)
					LNSD <- (log((theSD/theMean)^2 + 1) )^0.5
					lnormmatrix[i,(NumRep + 1):length(lnormmatrix[1,])] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
				} else {
					lnormmatrix[i,(NumRep + 1):length(lnormmatrix[1,])] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
					TBC <- sample(changing,1)
					if (TBC < 1) {
						Group1 <- c(Group1, i)
					}
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
	

	#Induce library size differences
	changes <- 2^(seq(-MaxLibSizelog2FC,0,0.001) )
	thesechanges <- sample(changes,length(lnormmatrix[1,]),replace=TRUE)
	for (i in seq_along(lnormmatrix[1,])) {
		lnormmatrix[,i] <- lnormmatrix[,i] * thesechanges[i]
	}
	lnormmatrix <-  as.matrix(floor(lnormmatrix))
	colnames(lnormmatrix) <- paste("Sample",1:ncol(lnormmatrix),sep="_")
	listing <- list(lnormmatrix,as.numeric(sort(tobechanged)), as.numeric(sort(Group1)), as.numeric(sort( tobechanged[!(tobechanged %in% Group1)])))
	results <- setNames(listing, c("data", "DiffList", "CorGroup1", "CorGroup2"))
	return (results)
}

NBSim <- function(thisdata_ori, NumRep=5, NumDiff = 2000, NumFea = 20000, showinfo=FALSE, MaxLibSizelog2FC=0.5, DEGlog2FC="Auto") {
	#RNA-seq data simulation for DEG analysis evaluation
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	#Capture distribution from the dataset for simulation
	options(warn=-1)
	#This step makes NBSim run slow, compared to the other 3 distributions, so, we will only model at most 10000 of the genes.
	Keep <- 1:nrow(thisdata_ori)
	if (nrow(thisdata_ori) > 10000) {
		Keep <- sample(Keep,10000)
	}
	
	MeanList <- vector(mode="numeric",length(Keep))
	d <- vector(mode="numeric",length(Keep))
	
	for (i in seq_along(thisdata_ori[Keep,1])){
		x <- as.numeric(thisdata_ori[i,])
		MeanList[i] <- mean(x)
		if (length(x) > 2) {
			d[i] <- as.numeric(unlist((glm.nb(x~1))[[24]]))
		} else {
			d[i] <- NA
		}
	}
	MeanList <- MeanList[!is.na(d)]
	d <- d[!is.na(d)]
	
	Fit <- LinearRegression(log(MeanList),log(d))
	FindDispersion <- function(inputmean) {
		return (exp(log(inputmean) * Fit$coefficients[[2]] + Fit$coefficients[[1]] ))
	}
	
	nbmatrix <- matrix(0, ncol=(2 * NumRep), nrow=NumFea)
	RN <- vector(mode="character", NumFea)
	for (i in 1:NumFea) {
		RN[i] <- paste("Gene",i,sep="")
	}
	rownames(nbmatrix) <- RN
	
	
	#Obtain Significant fold change threshold
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
	#FC of 5 is sure to be significant, so there is no need to search for boundary. If it ever happens to be not so, it also serves as safety for the program, since FC > 5  not being significant doesn't make sense.
	maxBound <- 5
	midBound <- round((minBound + maxBound)/2,4)
	themean <- median(MeanList)
	theDis <- FindDispersion(themean)
	
	NR2 <- NumRep
	design <- matrix(nrow=(NR2 * 2), ncol=2)
	colnames(design) <- c("SampleInfo", "Gam")
	design[,1] <- c(rep("Normal",NR2),rep("Tumor",NR2))
	design <- as.data.frame(design)
	pvstore <- vector(mode="numeric",100)
	if (DEGlog2FC != "Auto") {
		SigFC <- log(2^DEGlog2FC)
		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find its stddev.
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
	} else { 
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
	}
	#Non differential features will have mean of 0 and stddev of SigFC
	changes <- sort(rnorm(NumFea,0,FCstddev))
	begin <- round(NumDiff/2,0)
	ending <- round(NumFea - NumDiff/2,0)
	unchanged <- changes[begin:ending]
	changing <- c(changes[1:begin],changes[ending:length(changes)])

	Group1 <- c()
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
					if (TBC > 1) {
						Group1 <- c(Group1, i)
					}
					themean <- themean * TBC
					theDis <- FindDispersion(themean)
					nbmatrix[i,(NumRep + 1):length(nbmatrix[1,])] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
				} else {
					nbmatrix[i,(NumRep + 1):length(nbmatrix[1,])] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
					TBC <- sample(changing,1)
					if (TBC < 1) {
						Group1 <- c(Group1, i)
					}
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

	
	#Induce library size differences
	changes <- 2^(seq(-MaxLibSizelog2FC,0,0.001) )
	thesechanges <- sample(changes,length(nbmatrix[1,]),replace=TRUE)

	for (i in seq_along(nbmatrix[1,])) {
		nbmatrix[,i] <- nbmatrix[,i] * thesechanges[i]
	}

	nbmatrix <-  as.matrix(floor(nbmatrix))
	colnames(nbmatrix) <- paste("Sample",1:ncol(nbmatrix),sep="_")
	listing <- list(nbmatrix,as.numeric(sort(tobechanged)), as.numeric(sort(Group1)), as.numeric(sort( tobechanged[!(tobechanged %in% Group1)])))
	results <- setNames(listing, c("data", "DiffList", "CorGroup1", "CorGroup2"))

	return (results)
}
