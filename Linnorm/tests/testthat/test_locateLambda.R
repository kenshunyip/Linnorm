library(Linnorm)

context("Locate Lambda")
data(SEQC)
SEQC <- tXPM(SEQC)
SEQC <- SEQC[,colSums(SEQC != 0) == nrow(SEQC)]
SEQC <- SEQC[,order(colMeans(SEQC))]

LocateLambdaSEQC <- LocateLambda(SEQC, 100, 60000000)
LocateLambda_legacySEQC <- LocateLambda_legacy(SEQC, 100, 60000000)

PDSEQC <- abs((LocateLambdaSEQC - LocateLambda_legacySEQC)/mean(c(LocateLambda_legacySEQC,LocateLambdaSEQC))) * 100
AcceptableSEQC <- PDSEQC < 5

test_that("LocateLambda SEQC", {
	expect_true(AcceptableSEQC)
})

data(LIHC)
LIHC <- t(LIHC/1000000)
LIHC <- LIHC[,colSums(LIHC != 0) == nrow(LIHC)]
LIHC <- LIHC[,order(colMeans(LIHC))]

LocateLambdaLIHC <- LocateLambda(LIHC, 100, 60000000)
LocateLambda_legacyLIHC <- LocateLambda_legacy(LIHC, 100, 60000000)

PDLIHC <- abs((LocateLambdaLIHC - LocateLambda_legacyLIHC)/mean(c(LocateLambda_legacyLIHC,LocateLambdaLIHC))) * 100
AcceptableLIHC <- PDLIHC < 5

test_that("LocateLambda LIHC", {
	expect_true(AcceptableLIHC)
})
