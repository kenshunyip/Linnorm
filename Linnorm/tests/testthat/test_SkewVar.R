library(Linnorm)
library(moments)
library(matrixStats)

context("V and S lambda calculation check")

SkewVar_Raw <- function(GeneExp, lambda) {
	GeneExp <- log1p(GeneExp * lambda)
	RM <- as.numeric(rowMeans(GeneExp))
	RSD <- as.numeric(rowSDs(GeneExp))
	RSkew <- as.numeric(apply(GeneExp,1,skewness))
	
	SDLR <- LinearRegression(RM, RSD)
	SkewLR <- LinearRegression(RM, RSkew)
	# m <- results$coefficients[[2]]
	# c <- results$coefficients[[1]]
	SkewIntercept <- -SkewLR$coefficients[[1]]/SkewLR$coefficients[[2]]
	SkewIntegral <- 0
	if (SkewIntercept > min(RM) && SkewIntercept < max(RM)) {
		SkewIntegral <- (abs(((SkewLR$coefficients[[2]] * (SkewIntercept^2)/2) + (SkewLR$coefficients[[1]] * SkewIntercept)) - ((SkewLR$coefficients[[2]] * (min(RM)^2)/2) + (SkewLR$coefficients[[1]] * min(RM))) ) + abs(((SkewLR$coefficients[[2]] * (max(RM)^2)/2) + (SkewLR$coefficients[[1]] * max(RM))) - ((SkewLR$coefficients[[2]] * (SkewIntercept^2)/2) + (SkewLR$coefficients[[1]] * SkewIntercept)) ) )/(max(RM) - min(RM))
	} else {
		SkewIntegral <- abs(((SkewLR$coefficients[[2]] * (max(RM)^2)/2) + (SkewLR$coefficients[[1]] * max(RM))) - ((SkewLR$coefficients[[2]] * (min(RM)^2)/2) + (SkewLR$coefficients[[1]] * min(RM))) )/(max(RM) - min(RM))
	}
	return ((log1p(abs(SDLR$coefficients[[2]])) + 1)^2 + (log1p(SkewIntegral) + 1)^2)
}
SkewVar_Raw2 <- function(GeneExp, lambda) {
	GeneExp <- log1p(GeneExp * lambda)
	RM <- rowMeans(GeneExp)
	RSD <- rowSDs(GeneExp)
	RSkew <- apply(GeneExp,1,skewness)
	
	SDLR <- LinearRegression(RM, RSD)
	SkewLR <- LinearRegression(RM, RSkew)
	# m <- results$coefficients[[2]]
	# c <- results$coefficients[[1]]
	SkewIntercept <- -SkewLR$coefficients[[1]]/SkewLR$coefficients[[2]]
	SkewIntegral <- 0
	if (SkewIntercept > min(RM) && SkewIntercept < max(RM)) {
		SkewIntegral <- (abs((SkewIntercept - min(RM)) * ( (SkewLR$coefficients[[2]]/2) * (SkewIntercept + min(RM)) + SkewLR$coefficients[[1]])) + abs((max(RM) - SkewIntercept) * ( (SkewLR$coefficients[[2]]/2) * (max(RM) + SkewIntercept) + SkewLR$coefficients[[1]])))/(max(RM) - min(RM))
	} else {
		SkewIntegral <- abs((SkewLR$coefficients[[2]]/2) * (max(RM) + min(RM)) + SkewLR$coefficients[[1]])
	}
	return ((log1p(abs(SDLR$coefficients[[2]])) + 1)^2 + (log1p(SkewIntegral) + 1)^2)
}


data(LIHC)
LIHC <- LIHC/1000000
LIHC <- LIHC[rowSums(LIHC != 0) == ncol(LIHC), ]
LIHC <- LIHC[order(rowMeans(LIHC)),]

data(SEQC)
SEQC <- t(tXPM(SEQC))
SEQC <- SEQC[rowSums(SEQC != 0) == ncol(SEQC), ]
SEQC <- SEQC[order(rowMeans(SEQC)),]

LIHC2Raw2 <- SkewVar_Raw2(LIHC,1000000)
LIHC2 <- SkewVar(LIHC,1000000)
PDLIHC <- abs((LIHC2Raw2 - LIHC2)/mean(c(LIHC2,LIHC2Raw2))) * 100
AcceptableLIHC <- PDLIHC < 5
SEQC2Raw2 <- SkewVar_Raw2(SEQC,1000000)
SEQC2 <- SkewVar(SEQC,1000000)
PDSEQC <- abs((SEQC2Raw2 - SEQC2)/mean(c(SEQC2,SEQC2Raw2))) * 100
AcceptableSEQC <- PDSEQC < 5

test_that("SkewVar", {
	expect_true(AcceptableLIHC)
	expect_true(AcceptableSEQC)
})

LIHC2 <- SkewVar2(t(LIHC),1000000)
SEQC2 <- SkewVar2(t(SEQC),1000000)

PDLIHC <- abs((LIHC2Raw2 - LIHC2[1])/mean(c(LIHC2[1],LIHC2Raw2))) * 100
AcceptableLIHC <- PDLIHC < 5
PDSEQC <- abs((SEQC2Raw2 - SEQC2[1])/mean(c(SEQC2[1],SEQC2Raw2))) * 100
AcceptableSEQC <- PDSEQC < 5

test_that("SkewVar2", {
	expect_true(AcceptableLIHC)
	expect_true(AcceptableSEQC)
})

test_that("Integral simplification done correctly.", {
	expect_equal(SkewVar_Raw2(LIHC,1000000), SkewVar_Raw(LIHC,1000000))
	expect_equal(SkewVar_Raw2(SEQC,1000000), SkewVar_Raw(SEQC,1000000))
})

