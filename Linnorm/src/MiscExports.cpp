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

// [[Rcpp::export]]
double LocateLambda2( arma::mat GeneExp, double search_exponent);
RcppExport SEXP LocateLambda2(SEXP xSEXP,SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type GeneExp(xSEXP);
	Rcpp::traits::input_parameter< double >::type search_exponent(ySEXP);
    return Rcpp::wrap(LocateLambda2(GeneExp,search_exponent));
END_RCPP
}
