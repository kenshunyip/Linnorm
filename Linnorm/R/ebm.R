#Copyright (c) 2015, Institute for Systems Biology
#
#Permission is hereby granted, free of charge, to any person obtaining
#a copy of this software and associated documentation files (the
#"Software"), to deal in the Software without restriction, including
#without limitation the rights to use, copy, modify, merge, publish,
#distribute, sublicense, and/or sell copies of the Software, and to
#permit persons to whom the Software is furnished to do so, subject to
#the following conditions:

#The above copyright notice and this permission notice shall be
#included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Author: William Poole
# Ported to R: David L Gibbs
# Email: dgibbs@systemsbiology.org / william.poole@systemsbiology.org / tknijnen@systemsbiology.org
# Created: June 2015
#
# This file is edited for Linnorm.
# Edited: November 2016
# Editor: (Ken) Shun Hang Yip
# Editor Email: shunyip@bu.edu

popZScore <- function(data_vector) {
	.Call(popZScoreCPP, data_vector)
}

#Input: raw data vector (of one variable) with no missing samples. May be a list or an array.
#Output Transforemd data vector w.
transformData <- function(data_vector) {
    s = popZScore(data_vector)
    distr = ecdf(s)
    sapply(s, function(a) -2*log(distr(a)))
}

#Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array.
#       Note: Method does not deal with missing values within the data.
#Output: An m x m matrix of pairwise covariances between transformed raw data vectors
calculateCovariances <- function(data_matrix){
    transformed_data_matrix = apply(data_matrix, MARGIN=1, FUN=transformData)
    covar_matrix = cov(transformed_data_matrix)
    covar_matrix
}


#Input: A m x m numpy array of covariances between transformed data vectors and a vector of m p-values to combine.
#Output: A combined P-value.
#        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
combinePValues <- function(covar_matrix, p_values, extra_info = FALSE){
    N = ncol(covar_matrix) # number of samples
    df_fisher = 2.0*N
    Expected  = 2.0*N
    cov_sum <- (2*sum(covar_matrix[lower.tri(covar_matrix, diag=FALSE)]))
    Var = 4.0*N+cov_sum
    c = Var/(2.0*Expected)
    df_brown = (2.0*Expected^2)/Var
    if (df_brown > df_fisher) {
        df_brown = df_fisher
        c = 1.0
    }
	Answer <- rep(1, nrow(p_values))
	for (i in 1:nrow(p_values)) {
		x = 2.0*sum( -log(p_values[i,]) )
		Answer[i] <- pchisq(df=df_brown, q=x/c, lower.tail=FALSE)
	}
	return (Answer)
}

#Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array
#       A vector of m P-values to combine. May be a list or of type numpy.array.
#Output: A combined P-value.
#        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
empiricalBrownsMethod <- function(data_matrix, p_values, extra_info = FALSE) {
  # inputs must be numeric
    covar_matrix = calculateCovariances(data_matrix)
    return(combinePValues(covar_matrix, p_values, extra_info))
}
