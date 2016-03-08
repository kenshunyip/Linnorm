/*
The MIT License (MIT) 
Copyright (c) <2016> <Shun Hang Yip>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double gammaShape (arma::vec vec2);
RcppExport SEXP gammaShape(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type vec2(xSEXP);
    return Rcpp::wrap(gammaShape(vec2));
END_RCPP
}


// [[Rcpp::export]]
double LocateLambda(const arma::mat& GeneExp, double search_exponent);
RcppExport SEXP LocateLambda(SEXP xSEXP,SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type GeneExp(xSEXP);
	Rcpp::traits::input_parameter< double >::type search_exponent(ySEXP);
    return Rcpp::wrap(LocateLambda(GeneExp,search_exponent));
END_RCPP
}
