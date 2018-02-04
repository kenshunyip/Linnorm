library(Linnorm)
library(moments)
library(matrixStats)

context("Moment calculation check")

#Linnorm estimates sd and skewness using a one pass algorithm.
#So, we expect an average of less than 5% difference between true values and estimated values, instead of expecting them to be equal.
#Initialize datasets
data(LIHC)
LIHC <- LIHC/1000000
LIHC <- LIHC[rowSums(LIHC != 0) == ncol(LIHC), ]
LIHC <- t(LIHC[order(rowMeans(LIHC)),])
set.seed(123456)
LIHCTested <- sample(1:ncol(LIHC),10)

data(SEQC)
SEQC <- tXPM(as.matrix(SEQC))
SEQC <- SEQC[,colSums(SEQC != 0) == nrow(SEQC)]
SEQC <- SEQC[,order(colMeans(SEQC))]
set.seed(123456)
SEQCTested <- sample(1:ncol(SEQC),10)


#colSDs
AnswerLIHC <- colSds(LIHC)
AnswerSEQC <- colSds(SEQC)

LinnormAnswerLIHC <- colSDs(LIHC)
LinnormAnswerSEQC <- colSDs(SEQC)

AcceptableSD <- 0
for (i in 1:ncol(SEQC)) {
	AcceptableSD <- AcceptableSD + (abs((AnswerSEQC[i] - LinnormAnswerSEQC[i])/mean(c(LinnormAnswerSEQC[i],AnswerSEQC[i]))) * 100)
}
Acceptable <- AcceptableSD/ncol(SEQC) < 5
test_that("colSDs SEQC mean test", {
	expect_true(Acceptable)
})


AcceptableSD <- 0
for (i in 1:ncol(LIHC)) {
	AcceptableSD <- AcceptableSD + (abs((AnswerLIHC[i] - LinnormAnswerLIHC[i])/mean(c(LinnormAnswerLIHC[i],AnswerLIHC[i]))) * 100)
}
Acceptable <- AcceptableSD/ncol(LIHC) < 5
test_that("colSDs LIHC mean test", {
	expect_true(Acceptable)
})

#NZcolMeans
NZcolMeans2 <- function(x) {
	answer <- rep(0, ncol(x))
	for (i in 1:ncol(x)) {
		answer[i] <- mean(x[x[,i] != 0,i])
	}
	return(answer)
}

AnswerLIHC <- NZcolMeans2(LIHC)
AnswerSEQC <- NZcolMeans2(SEQC)
LinnormAnswerLIHC <- NZcolMeans(LIHC)
LinnormAnswerSEQC <- NZcolMeans(SEQC)

AcceptableMean <- 0
for (i in 1:ncol(SEQC)) {
	AcceptableMean <- AcceptableMean + (abs((AnswerSEQC[i] - LinnormAnswerSEQC[i])/mean(c(LinnormAnswerSEQC[i],AnswerSEQC[i]))) * 100)
}
Acceptable <- AcceptableMean/ncol(SEQC) < 5
test_that("NZcolMeans SEQC mean test", {
	expect_true(Acceptable)
})

AcceptableMean <- 0
for (i in 1:ncol(LIHC)) {
	AcceptableMean <- AcceptableMean + (abs((AnswerLIHC[i] - LinnormAnswerLIHC[i])/mean(c(LinnormAnswerLIHC[i],AnswerLIHC[i]))) * 100)
}
Acceptable <- AcceptableMean/ncol(LIHC) < 5
test_that("NZcolMeans LIHC mean test", {
	expect_true(Acceptable)
})



#NZcolMeanSD
NZcolMeanSD2 <- function(x) {
	answer <- matrix(nrow=2, ncol=ncol(x))
	for (i in 1:ncol(x)) {
		answer[1,i] <- mean(x[x[,i] != 0,i])
		answer[2,i] <- sd(x[x[,i] != 0,i])
	}
	return(answer)
}

AnswerLIHC <- NZcolMeanSD2(LIHC)
AnswerSEQC <- NZcolMeanSD2(SEQC)
LinnormAnswerLIHC <- NZcolMeanSD(LIHC)
LinnormAnswerSEQC <- NZcolMeanSD(SEQC)

AcceptableMean <- 0
AcceptableSD <- 0
for (i in 1:ncol(SEQC)) {
	AcceptableMean <- AcceptableMean + (abs((AnswerSEQC[1,i] - LinnormAnswerSEQC[1,i])/mean(c(LinnormAnswerSEQC[1,i],AnswerSEQC[1,i]))) * 100)
	AcceptableSD <- AcceptableSD + (abs((AnswerSEQC[2,i] - LinnormAnswerSEQC[2,i])/mean(c(LinnormAnswerSEQC[2,i],AnswerSEQC[2,i]))) * 100)
}
Acceptable <- AcceptableMean/ncol(SEQC) < 5
test_that("NZcolMeanSD SEQC mean test", {
	expect_true(Acceptable)
})
Acceptable <- AcceptableSD/ncol(SEQC) < 5
test_that("NZcolMeanSD SEQC sd test", {
	expect_true(Acceptable)
})


AcceptableMean <- 0
AcceptableSD <- 0
for (i in 1:ncol(LIHC)) {
	AcceptableMean <- AcceptableMean + (abs((AnswerLIHC[1,i] - LinnormAnswerLIHC[1,i])/mean(c(LinnormAnswerLIHC[1,i],AnswerLIHC[1,i]))) * 100)
	AcceptableSD <- AcceptableSD + (abs((AnswerLIHC[2,i] - LinnormAnswerLIHC[2,i])/mean(c(LinnormAnswerLIHC[2,i],AnswerLIHC[2,i]))) * 100)
}
Acceptable <- AcceptableMean/ncol(LIHC) < 5
test_that("NZcolMeanSD LIHC mean test", {
	expect_true(Acceptable)
})
Acceptable <- AcceptableSD/ncol(LIHC) < 5
test_that("NZcolMeanSD LIHC sd test", {
	expect_true(Acceptable)
})


#NZcolLogMeanSDSkew
NZcolLogMeanSDSkew2 <- function(x) {
	answer <- matrix(nrow=3, ncol=ncol(x))
	for (i in 1:ncol(x)) {
		thisdata <- log(x[x[,i] != 0,i])
		answer[1,i] <- mean(thisdata)
		answer[2,i] <- sd(thisdata)
		answer[3,i] <- skewness(thisdata)
	}
	return(answer)
}

AnswerLIHC <- NZcolLogMeanSDSkew2(LIHC)
AnswerSEQC <- NZcolLogMeanSDSkew2(SEQC)
LinnormAnswerLIHC <- NZcolLogMeanSDSkew(LIHC)
LinnormAnswerSEQC <- NZcolLogMeanSDSkew(SEQC)

AcceptableMean <- 0
AcceptableSD <- 0
AcceptableSkew <- 0
for (i in 1:ncol(SEQC)) {
	AcceptableMean <- AcceptableMean + (abs((AnswerSEQC[1,i] - LinnormAnswerSEQC[1,i])/mean(c(LinnormAnswerSEQC[1,i],AnswerSEQC[1,i]))) * 100)
	AcceptableSD <- AcceptableSD + (abs((AnswerSEQC[2,i] - LinnormAnswerSEQC[2,i])/mean(c(LinnormAnswerSEQC[2,i],AnswerSEQC[2,i]))) * 100)
	AcceptableSkew <- AcceptableSkew + (abs((AnswerSEQC[3,i] - LinnormAnswerSEQC[3,i])/mean(c(LinnormAnswerSEQC[3,i],AnswerSEQC[3,i]))) * 100)
}
Acceptable <- AcceptableMean/ncol(SEQC) < 5
test_that("NZcolLogMeanSDSkew SEQC mean test", {
	expect_true(Acceptable)
})
Acceptable <- AcceptableSD/ncol(SEQC) < 5
test_that("NZcolLogMeanSDSkew SEQC sd test", {
	expect_true(Acceptable)
})
Acceptable <- AcceptableSkew/ncol(SEQC) < 5
test_that("NZcolLogMeanSDSkew SEQC skew test", {
	expect_true(Acceptable)
})

AcceptableMean <- 0
AcceptableSD <- 0
AcceptableSkew <- 0
for (i in 1:ncol(LIHC)) {
	AcceptableMean <- AcceptableMean + (abs((AnswerLIHC[1,i] - LinnormAnswerLIHC[1,i])/mean(c(LinnormAnswerLIHC[1,i],AnswerLIHC[1,i]))) * 100)
	AcceptableSD <- AcceptableSD + (abs((AnswerLIHC[2,i] - LinnormAnswerLIHC[2,i])/mean(c(LinnormAnswerLIHC[2,i],AnswerLIHC[2,i]))) * 100)
	AcceptableSkew <- AcceptableSkew + (abs((AnswerLIHC[3,i] - LinnormAnswerLIHC[3,i])/mean(c(LinnormAnswerLIHC[3,i],AnswerLIHC[3,i]))) * 100)
}
Acceptable <- AcceptableMean/ncol(LIHC) < 5
test_that("NZcolLogMeanSDSkew LIHC mean test", {
	expect_true(Acceptable)
})
Acceptable <- AcceptableSD/ncol(LIHC) < 5
test_that("NZcolLogMeanSDSkew LIHC sd test", {
	expect_true(Acceptable)
})
Acceptable <- AcceptableSkew/ncol(LIHC) < 5
test_that("NZcolLogMeanSDSkew LIHC skew test", {
	expect_true(Acceptable)
})
