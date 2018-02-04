library(Linnorm)

context("Linear Regression check")

#Initialize datasets
a <- rnorm(100)
b <- rnorm(100)

LR <- LinearRegression(a,b)
LR2 <- lm(b~a)

d <- 1:100
e <- 110:11
LR3 <- LinearRegression(d,e)
LR4 <- lm(e~d)
test_that("Linear Regression", {
	expect_equal(LR$coefficients[[2]], LR$coefficients[[2]])
	expect_equal(LR$coefficients[[1]], LR$coefficients[[1]])
	expect_equal(LR3$coefficients[[2]], LR4$coefficients[[2]])
	expect_equal(LR3$coefficients[[1]], LR4$coefficients[[1]])
	
})
