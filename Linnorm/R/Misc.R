#Convert Dataset into XPM
XPM <- function(x) {
    .Call(XPMCpp, x)
}
#Convert Dataset into XPM and transpose
tXPM <- function(x) {
    .Call(tXPMCpp, x)
}

#Find low count gene filtering threshold
FindLCT <- function(datamatrix, Multy,showinfo) {
	#Find low count gene filtering threshold
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	#Obtain mean and stdev for log(Relative Expression * estimated total count + 1) dataset
	MeanSD <- NZcolLog1pMeanSD(datamatrix,Multy)
	LC_Threshold <- 0
	MeanSD <- MeanSD[,which(!is.nan(MeanSD[2,]))]
	#Lowest expressing 1/3 of the genes will be used for the calculation of slope
	Portion <- 1/3
	MeanOrder <- order(MeanSD[1,])
	Slope <- 1
	#Loop until Slope is negative for 3 times in a roll or thersohld is too big
	numNegative <- 0
	while(numNegative < 3 && LC_Threshold < 0.98) {
		LC_Threshold <- LC_Threshold + 0.01
		Range <- floor(ncol(MeanSD) * LC_Threshold + 1):ncol(MeanSD)
		Range <- Range[1:floor(length(Range) * Portion + 1)]
		Slope <- getSlope(MeanSD[1,MeanOrder[Range]],MeanSD[2,MeanOrder[Range]])
		if (Slope < 0) {
			numNegative <- numNegative + 1
		} else {
			numNegative <- 0
		}
	}
	LC_Threshold <- LC_Threshold - ((numNegative - 1)/100)
	return (round(LC_Threshold,2))
}
FindLCT_DI <- function(datamatrix, showinfo) {
	#Find low count gene filtering threshold for Data imputation funciton, where the dataset is Linnorm transformed.
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	#Obtain mean and stdev
	MeanSD <- colMeanSD(datamatrix)
	LC_Threshold <- 0
	MeanSD <- MeanSD[,which(!is.nan(MeanSD[2,]))]
	MeanOrder <- order(MeanSD[1,])
	Slope <- 1
	while(Slope > 0 && LC_Threshold < 1) {
		LC_Threshold <- LC_Threshold + 0.01
		Range <- floor(ncol(MeanSD) * LC_Threshold + 1):ncol(MeanSD)
		Slope <- getSlope(MeanSD[1,MeanOrder[Range]],MeanSD[2,MeanOrder[Range]])
	}
	return (round(LC_Threshold,2))
}

#Filter dataset and obtain model genes
FirstFilter <- function(x, minNonZeroPortion, L_F_p = 0.25, L_F_LC_Genes = 0.01, L_F_HC_Genes = 0.01, spikein = NULL) {
	#Filter dataset and obtain model genes for calculation of lambda
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	#Obtain mean, stdev and skewness
	MeanSDSkew <- NZcolLogMeanSDSkew(x)
	Keep <- which(!is.nan(MeanSDSkew[3,]))
	
	#Sort data and filter with L_F_LC_Genes and minNonZeroPortion
	if (minNonZeroPortion == 0 || minNonZeroPortion == 1) {
		Keep <- Keep[which(colSums(x[,Keep] != 0) >= nrow(x) * minNonZeroPortion)]
	} else {
		Keep <- Keep[which(colSums(x[,Keep] != 0) > nrow(x) * minNonZeroPortion)]
	}
	
	Keep <- Keep[order(MeanSDSkew[1,Keep], decreasing = FALSE)]
	
	Start <- floor(length(Keep) * L_F_LC_Genes + 1)
	End <- length(Keep) - floor(length(Keep) * L_F_HC_Genes)
	Keep <- Keep[Start:End]
	
	x <- x[,Keep]
	MeanSDSkew <- MeanSDSkew[,Keep]
	
	#Check the number if spike in genes and whether they are provided
	spikeino <- spikein
	spikein <- spikein[which(spikein %in% colnames(x))]
	if (length(spikein) < 3 && length(spikeino) != 0) {
		spikein = NULL
		warning("Too many Spikein are filtered. They will not be utilized.")
	}
	
	#Initialize object for storing genes that will be retained
	allStableGenes <- 0
	
	#LOWESS fit, Use precision weight as weight
	logitit <- loessFit(MeanSDSkew[2,],MeanSDSkew[1,], weights=1/MeanSDSkew[2,]^2)
	Keep <- which(logitit$fitted > 0)
	LogFit <- logitit$fitted[Keep]
	x <- x[,Keep]
	MeanSDSkew <- MeanSDSkew[,Keep]
	SDRatio <- log(as.numeric(MeanSDSkew[2,]/LogFit))
	
	#normalize SDRatio
	LR <- LinearRegression(MeanSDSkew[1,],SDRatio)
	Residual <- SDRatio - LR$coefficients[[2]] * MeanSDSkew[1,] - LR$coefficients[[1]]
	LR2 <- LinearRegression(MeanSDSkew[1,],abs(Residual))
	SDRatio <- SDRatio * (LR2$coefficients[[2]] * MeanSDSkew[1,1] + LR2$coefficients[[1]])/(LR2$coefficients[[2]] * MeanSDSkew[1,] + LR2$coefficients[[1]])
	
	#LOWESS fit for skewness
	SkewResidual <- loessFit(MeanSDSkew[3,],MeanSDSkew[1,])
	SkewResidual <- SkewResidual$residuals
	
	#Object for storing pvalues from stdev and skewness. Column 1: Stdev p values. Column 2: Skew p values.
	pvalueMatrix <- matrix(nrow=ncol(MeanSDSkew), ncol=2)
	
	
	if (length(spikein) < 3) {
	#degree of freedom is 2 (using loess() and predict() to estimate df is very slow. Here, we will just assume it to be n-2 to save time. Given hundreds to tens of thousands of features in RNA-seq dataset, this df should be a good enough estimate.)
		SDnoOutlier <- SDRatio[!SDRatio %in% boxplot.stats(SDRatio)$out]
		TheMean <- mean(SDnoOutlier)
		tdeno <- sqrt(sum((SDnoOutlier - TheMean)^2)/(length(SDnoOutlier) - 2))
		pvalueMatrix[,1] <- 2 * pt(abs((SDRatio - TheMean)/tdeno), df = length(SDnoOutlier) - 2, lower.tail = FALSE)
		
		SkewnoOutlier <- SkewResidual[!SkewResidual %in% boxplot.stats(SkewResidual)$out]
		TheMean <- mean(SkewnoOutlier)
		tdeno <- sqrt(sum((SkewnoOutlier - TheMean)^2)/(length(SkewnoOutlier) - 2))
		pvalueMatrix[,2] <- 2 * pt(abs((SkewResidual - TheMean)/tdeno),df = length(SkewnoOutlier) - 2, lower.tail = FALSE)
	} else {
	#degree of freedom is 2 (using loess() and predict() to estimate df is very slow. Here, we will just assume it to be n-2 to save time. Given hundreds to tens of thousands of features in RNA-seq dataset, this df should be a good enough estimate.)
		spikes <- which(colnames(x) %in% spikein)
		SDnoOutlier <- SDRatio[spikes]
		TheMean <- mean(SDnoOutlier)
		tdeno <- sqrt(sum((SDnoOutlier - TheMean)^2)/(length(SDnoOutlier) - 2))
		pvalueMatrix[,1] <- 2 * pt(abs((SDRatio - TheMean)/tdeno), df = length(SDnoOutlier) - 2, lower.tail = FALSE)
		SkewnoOutlier <- SkewResidual[spikes]
		TheMean <- mean(SkewnoOutlier)
		tdeno <- sqrt(sum((SkewnoOutlier - TheMean)^2)/(length(SkewnoOutlier) - 2))
		pvalueMatrix[,2] <- 2 * pt(abs((SkewResidual - TheMean)/tdeno),df = length(SkewnoOutlier) - 2, lower.tail = FALSE)
	}
	
	#Combine p values
	#combinedPvalues <- apply(pvalueMatrix,1,FisherMethod)
	combinedPvalues <- empiricalBrownsMethod(MeanSDSkew[2:3,],pvalueMatrix)
	
	#Obtain stable genes
	allStableGenes <- which(combinedPvalues > L_F_p)
	
	#Safety, need 100 genes at least
	while (length(allStableGenes) < 100) {
		porder <- order(combinedPvalues, decreasing=TRUE)
		allStableGenes <- porder[1:100]
	}
	#Output
	return (x[,allStableGenes])
}

BatchEffectLinnorm1 <- function(x, minNonZeroPortion, BE_F_LC_Genes = 0.25,BE_F_HC_Genes = 0.05, BE_F_p = 0.5, BE_strength = 0.25, spikein = NULL) {
	#Normalization for batch effect.
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	
	#Save a copy of the raw matrix. x2 will not be filtered and used for normalization. x will be filtered and used as the model.
	x2 <- x
	
#Filter dataset and obtain model genes
	#Obtain mean, stdev and skewness
	MeanSDSkew <- NZcolLogMeanSDSkew(x)
	
	Keep <- which(!is.nan(MeanSDSkew[3,]))
	Keep <- Keep[order(MeanSDSkew[1,Keep], decreasing = FALSE)]
	
	#Sort data and filter with BE_F_LC_Genes and minNonZeroPortion
	if (minNonZeroPortion == 0 || minNonZeroPortion == 1) {
		Keep <- Keep[which(colSums(x[,Keep] != 0) >= nrow(x) * minNonZeroPortion)]
	} else {
		Keep <- Keep[which(colSums(x[,Keep] != 0) > nrow(x) * minNonZeroPortion)]
	}
	
	Start <- floor(length(Keep) * BE_F_LC_Genes + 1)
	End <- length(Keep) - floor(length(Keep) * BE_F_HC_Genes)
	Keep <- Keep[Start:End]
	
	x <- x[,Keep]
	MeanSDSkew <- MeanSDSkew[,Keep]
	
	#Check the number if spike in genes and whether they are provided
	spikeino <- spikein
	spikein <- spikein[which(spikein %in% colnames(x))]
	if (length(spikein) < 3 && length(spikeino) != 0) {
		spikein = NULL
		warning("Too many Spikein are filtered. They will not be utilized.")
	}
	
	#Initialize object for storing genes that will be retained
	allStableGenes <- 0
	
	#LOWESS fit, Use precision weight as weight
	logitit <- loessFit(MeanSDSkew[2,],MeanSDSkew[1,], weights=1/MeanSDSkew[2,]^2)	
	Keep <- which(logitit$fitted > 0)
	LogFit <- logitit$fitted[Keep]
	x <- x[,Keep]
	MeanSDSkew <- MeanSDSkew[,Keep]	
	SDRatio <- log(as.numeric(MeanSDSkew[2,]/LogFit))
	
	#normalize SDRatio
	LR <- LinearRegression(MeanSDSkew[1,],SDRatio)
	Residual <- SDRatio - (LR$coefficients[[2]] * MeanSDSkew[1,] + LR$coefficients[[1]])
	LR2 <- LinearRegression(MeanSDSkew[1,],abs(Residual))
	SDRatio <- SDRatio * (LR2$coefficients[[2]] * MeanSDSkew[1,1] + LR2$coefficients[[1]])/(LR2$coefficients[[2]] * MeanSDSkew[1,] + LR2$coefficients[[1]])
	
	#LOWESS fit for skewness
	SkewResidual <- loessFit(MeanSDSkew[3,],MeanSDSkew[1,])
	SkewResidual <- SkewResidual$residuals
	
	#Object for storing pvalues from stdev and skewness. Column 1: Stdev p values. Column 2: Skew p values.
	pvalueMatrix <- matrix(nrow=ncol(MeanSDSkew), ncol=2)
	
	if (length(spikein) < 3) {
	#degree of freedom is 2 (using loess and predict to estimate df is very slow. Here, we will just assume it to be n-2 to save time. Given hundreds to tens of thousands of features in RNA-seq dataset, this df should be a good enough estimate.)
		SDnoOutlier <- SDRatio[!SDRatio %in% boxplot.stats(SDRatio)$out]
		TheMean <- mean(SDnoOutlier)
		tdeno <- sqrt(sum((SDnoOutlier - TheMean)^2)/(length(SDnoOutlier) - 2))
		pvalueMatrix[,1] <- 2 * pt(abs((SDRatio - TheMean)/tdeno), df = length(SDnoOutlier) - 2, lower.tail = FALSE)
		
		SkewnoOutlier <- SkewResidual[!SkewResidual %in% boxplot.stats(SkewResidual)$out]
		TheMean <- mean(SkewnoOutlier)
		tdeno <- sqrt(sum((SkewnoOutlier - TheMean)^2)/(length(SkewnoOutlier) - 2))
		pvalueMatrix[,2] <- 2 * pt(abs((SkewResidual - TheMean)/tdeno),df = length(SkewnoOutlier) - 2, lower.tail = FALSE)
	} else {
	#degree of freedom is 2 (using loess() and predict() to estimate df is very slow. Here, we will just assume it to be n-2 to save time. Given hundreds to tens of thousands of features in RNA-seq dataset, this df should be a good enough estimate.)
		spikes <- which(colnames(x) %in% spikein)
		SDnoOutlier <- SDRatio[spikes]
		TheMean <- mean(SDnoOutlier)
		tdeno <- sqrt(sum((SDnoOutlier - TheMean)^2)/(length(SDnoOutlier) - 2))
		
		pvalueMatrix[,1] <- 2 * pt(abs((SDRatio - TheMean)/tdeno), df = length(SDnoOutlier) - 2, lower.tail = FALSE)
		
		
		SkewnoOutlier <- SkewResidual[spikes]
		TheMean <- mean(SkewnoOutlier)
		tdeno <- sqrt(sum((SkewnoOutlier - TheMean)^2)/(length(SkewnoOutlier) - 2))
		pvalueMatrix[,2] <- 2 * pt(abs((SkewResidual - TheMean)/tdeno),df = length(SkewnoOutlier) - 2, lower.tail = FALSE)
	}
	#Combine p values
	#combinedPvalues <- apply(pvalueMatrix,1,FisherMethod)
	combinedPvalues <- empiricalBrownsMethod(MeanSDSkew[2:3,],pvalueMatrix)
	combinedPvalues[is.na(combinedPvalues)] <- 0

	#Obtain stable genes
	allStableGenes <- which(combinedPvalues > BE_F_p)

	#Safety, need 100 genes at least
	while (length(allStableGenes) < 100) {
		porder <- order(combinedPvalues, decreasing=TRUE)
		allStableGenes <- porder[1:100]
	}
	CN <- colnames(x2)
	RN <- rownames(x2)
	
	#Using the model stable genes, perform normalization
	x2 <- BatchEffect2(x[,allStableGenes], x2, MeanSDSkew[1,allStableGenes], BE_strength)
	colnames(x2) <- CN
	rownames(x2) <- RN
	
	#Output
	return (x2)
}
BatchEffect2 <- function(x,y,z,z2) {
	.Call(BatchEffectCpp, x,y,z,z2)
}

#Linnorm's main funciton. Find optimal Lambda. This is implemented in C++.
LocateLambda <- function(x,y,z) {
    .Call(LocateLambdaCpp, x,y,z)
}
LocateLambda_legacy <- function(x,y,z) {
    .Call(LocateLambdaCpp_legacy, t(x),y,z)
}

#Legacy functions used for testing.
SkewVar <- function(x,y) {
    .Call(SkewVarCpp, x,y)
}
SkewVar2 <- function(x,y) {
    .Call(SkewVar2Cpp, x,y)
}
SkewAVar <- function(x,y) {
    .Call(SkewAVarCpp, x,y)
}

createUpperIndex <- function(colLength,TotalLength) {
	#Create index for parsing correlation matrix on the "upper triangle" vector.
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	index <- matrix(0,nrow=TotalLength,ncol=2)
	TL <- 1
	for (i in 2:colLength) {
		newmatrix <- matrix(i, ncol=2, nrow=(i-1))
		newmatrix[,1] <- seq(1,(i-1),1)
		index[TL:(TL+i-2),] <- newmatrix
		TL <- TL + i - 1
	}
	return (index)
}
createLowerIndex <- function(rowLength,TotalLength) {
	#Same as above, but for "lower triangle"
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	index <- matrix(0,nrow=TotalLength,ncol=2)
	TL <- 1
	for (i in 2:rowLength) {
		newmatrix <- matrix(i, ncol=2, nrow=(i-1))
		newmatrix[,2] <- seq(1,(i-1),1)
		index[TL:(TL+i-2),] <- newmatrix
		TL <- TL + i - 1
	}
	return (index)
}
UpperToMatrix <- function(datavalues,UpperIndex) {
	#Convert  "upper triangle" vector to matrix.
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	theMatrix <- matrix(1, ncol=max(UpperIndex[,2]), nrow=max(UpperIndex[,2]))
	upperrow <- order(UpperIndex[,1])
	lowerrow <- order(UpperIndex[,2])
	upperrowi <- 1
	upperrowj <- 1
	lowerrowi <- 1
	lowerrowj <- 1
	for (i in 1:ncol(theMatrix)) {
		if (upperrowj < length(upperrow)) {
			while (UpperIndex[upperrow[upperrowj],1] == i) {
				upperrowj <- upperrowj + 1
				if (upperrowj == length(upperrow)) {
					upperrowj <- upperrowj + 1
					break
				}
			}
			theMatrix[i,UpperIndex[upperrow[upperrowi:(upperrowj-1)],2]] <- datavalues[upperrow[upperrowi:(upperrowj-1)]]
		}
		if (lowerrowj < length(lowerrow)) {
			while (UpperIndex[lowerrow[lowerrowj],2] == i) {
				lowerrowj <- lowerrowj + 1
				if (lowerrowj == length(lowerrow)) {
					lowerrowj <- lowerrowj + 1
					break
				}
			}
			theMatrix[i,UpperIndex[lowerrow[lowerrowi:(lowerrowj-1)],1]] <- datavalues[lowerrow[lowerrowi:(lowerrowj-1)]]
		}
		upperrowi <- upperrowj
		lowerrowi <- lowerrowj
	}
	return (theMatrix)
}

areColors <- function(x) {
	#Check if input are colors
	sapply(x, function(X) {
		tryCatch(is.matrix(col2rgb(X)), 
		error = function(e) FALSE)
	})
}
