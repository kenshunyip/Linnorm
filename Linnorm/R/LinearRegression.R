#' One Pass Linear Regression.
#'
#' This function performs Linear Regression on the input data with a one pass algorithm implemented in C++. This is for users who only need m and c from the y=mx + c equation. Compared to the lm function, this function is much faster. 
#' @param x	Numeric vector. x values.
#' @param y	Numeric vector. corresponding y values.
#' @details  This function calculates m and c in the linear equestion, y = mx + c. 
#' @return This function returns a list with one object, "coefficients". The first element in this object is c; the second element is m in the y = mx + c equation.
#' @keywords Linear Regression slope constant
#' @export
#' @examples
#' x <- 1:10
#' y <- 1:10
#' results <- LinearRegression(x,y)
#' m <- results$coefficients[[2]]
#' c <- results$coefficients[[1]]
#' @import
#' Rcpp
#' RcppArmadillo
LinearRegression <- function(x,y) {
	#Linear regression function
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	if (length(x) != length(y)) {
		stop("x and y do not have the same length.")
	}
	#Call C++
	Answer <- LR(x,y)
	#Put results into list for output
	listing <- list(Answer[2], Answer[1])
	listing2 <- list(listing)
	result <- setNames(listing2, c("coefficients"))
	return (result)
}

LR <- function(x,y) {
	.Call(LinearRegressionCpp, x,y)
}

LinearRegressionZero <- function(x,y) {
	#Linear regression through the origin.
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	if (length(x) != length(y)) {
		stop("x and y do not have the same length.")
	}
	#Call C++
	Answer <- LRZ(x,y)
	#Put results into list for output
	listing <- list(Answer[2], Answer[1])
	listing2 <- list(listing)
	result <- setNames(listing2, c("coefficients"))
	return (result)
}

LRZ <- function(x,y) {
	.Call(LinearRegressionZeroCpp, x,y)
}

#' One Pass Linear Regression with fixed point.
#'
#' This function performs Linear Regression on the input data with a fixed point. It uses a one pass algorithm implemented in C++. This is for users who only need m and c from the y=mx + c equation. Compared to the lm function, this function is much faster. 
#' @param x	Numeric vector. x values.
#' @param y	Numeric vector. corresponding y values.
#' @param x1	Numeric. x coordinate of the fixed point.
#' @param y1	Numeric. y coordinate of the fixed point.
#' @details  This function calculates m and c in the linear equestion, y = mx + c. 
#' @return This function returns a list with one object, "coefficients". The first element in this object is c; the second element is m in the y = mx + c equation.
#' @keywords Linear Regression slope constant
#' @export
#' @examples
#' x <- 1:10
#' y <- 1:10
#' x1 <- 1
#' y1 <- 2
#' results <- LinearRegressionFP(x,y, x1, y1)
#' m <- results$coefficients[[2]]
#' c <- results$coefficients[[1]]
#' @import
#' Rcpp
#' RcppArmadillo
LinearRegressionFP <- function(x,y, x1,y1) {
	#Linear regression through a fixed point.
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	if (length(x) != length(y)) {
		stop("x and y do not have the same length.")
	}
	#Call C++
	Answer <- 0
	if (x1 == 0 && y1 == 0) {
		Answer <- LRZ(x,y)
	} else {
		Answer <- LRFP(x,y,x1,y1)
	}
	#Put results into list for output
	listing <- list(Answer[2], Answer[1])
	listing2 <- list(listing)
	result <- setNames(listing2, c("coefficients"))
	return (result)
}

LRFP <- function(x,y, x1,y1) {
	.Call(LinearRegressionFPCpp, x,y, x1,y1)
}

#Get Slope from x and y vectors
getSlope <- function(x,y) {
	.Call(getSlopeCpp, x,y)
}
