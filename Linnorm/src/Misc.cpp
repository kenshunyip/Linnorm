/*
The MIT License (MIT) 
Copyright (c) <2017> <(Ken) Shun Hang Yip> <shunyip@bu.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/
#include <R_ext/Rdynload.h>
#include <R.h>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <stdint.h>
// [[Rcpp::plugins(cpp11)]]


using namespace std;
using namespace Rcpp;

//Convert dataset into Relative Expression. Please note that even though they are called "Per million", they ceased to be per million during development.
SEXP XPMCpp (SEXP xSEXP){
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	arma::colvec RowSums = sum(GeneExp, 1);
	int_fast32_t i=0, n=0;
	arma::mat::iterator it_end = GeneExp.end();
	//One pass linear regression with one pass variance, skewness
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		*it = *it/RowSums.at(n);
		n++;
		if (n == int(GeneExp.n_rows)) {
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(GeneExp);
}
SEXP tXPMCpp (SEXP xSEXP){
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	arma::rowvec ColumnSums = sum(GeneExp, 0);
	int_fast32_t i=0, n=0;
	arma::mat::iterator it_end = GeneExp.end();
	//One pass linear regression with one pass variance, skewness
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		*it = *it/ColumnSums.at(i);
		n++;
		if (n == int(GeneExp.n_rows)) {
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(trans(GeneExp));
}

//The following three functions calculate the F(lambda) for the dataset given a lambda.
//This function is the original one written, and is used for unit test (together with the other two R functions.)
static double SkewVar (const arma::mat& GeneExp, const double& lambda2) {
	//vectors to store skewness, standard deviation and mean of each gene/feature.

	double thisSD, thisSkew;
	vector< double > meanvec;
	//Objects to store the sums of variables that will be needed to perform linear regression.
	double SumSkew = 0;
	double SumSD = 0;
	double SumMean = 0;
	double SumSkewMean = 0;
	double SumSDMean = 0;
	double SumMeanSq = 0;
	
	double mean = 0;
	double M2 = 0;
	double M3 = 0;
	double numData = 0;
	double delta, delta_n, term1;
		
	//One pass linear regression with one pass variance, skewness
	for (int_fast32_t i = 0; i < int_fast32_t(GeneExp.n_rows); i ++) {
		mean = 0;
		M2 = 0;
		M3 = 0;
		numData = 0;
		for (int_fast32_t n = 0; n < int_fast32_t(GeneExp.n_cols); n++) {
			if (GeneExp.at(i,n) != 0) {
				delta = log1p(GeneExp.at(i,n) * lambda2) - mean;
				delta_n = delta / (numData + 1);
				term1 = delta * delta_n * numData;
				mean = mean + delta_n;
				M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
				M2 = M2 + term1;
				numData++;
			}
		}
		//Here, calculate skewness and SD using 3rd and 2nd moments.
		thisSkew =  (sqrt(numData) * M3) / pow(M2,1.5);
		thisSD = sqrt(M2/(numData - 1));
		
		//Linear regression
		meanvec.push_back(mean);
		SumSkew += thisSkew;
		SumSD += thisSD;
		SumMean += meanvec.at(i);
		SumSkewMean += thisSkew * mean;
		SumSDMean += thisSD * mean;
		SumMeanSq += pow(meanvec.at(i),2);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	double skewM = (SumSkewMean * GeneExp.n_rows - SumSkew * SumMean)/(SumMeanSq * GeneExp.n_rows - pow(SumMean,2) ) ;
	double skewC = (SumSkew - skewM * SumMean)/GeneExp.n_rows;
	double SDevM = (SumSDMean * GeneExp.n_rows - SumSD * SumMean)/(SumMeanSq * GeneExp.n_rows - pow(SumMean,2) );
	
	//integral of the linear equation of skewness
	double Skewzerointercept = -skewC/skewM;
	double Skewintegral;
	skewM = skewM/2;
	if (Skewzerointercept > meanvec.at(0) && Skewzerointercept < meanvec.at(GeneExp.n_rows-1)) {
		Skewintegral = (abs(skewM * (meanvec.at(GeneExp.n_rows-1) + Skewzerointercept) + skewC) * (meanvec.at(GeneExp.n_rows-1) - Skewzerointercept) + abs(skewM * (meanvec.at(0) + Skewzerointercept) + skewC)* (Skewzerointercept - meanvec.at(0)))/(meanvec.at(GeneExp.n_rows-1) - meanvec.at(0));
	} else {
		Skewintegral = abs(skewM * (meanvec.at(GeneExp.n_rows-1) + meanvec.at(0)) + skewC) ;
	}
	return pow(log1p(abs(SDevM))+1,2) + pow(log1p(Skewintegral)+1,2);
}
SEXP SkewVarCpp(SEXP xSEXP,SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double lambda2 = Rcpp::as<double>(ySEXP);
	return Rcpp::wrap(SkewVar(GeneExp, lambda2));
}
//Optimized for speed
//This function is the same as SkewVar.
static double SkewVar1 (arma::mat& GeneExp, const double& lambda2) {
	//double makes the program slower, but it is needed for accurate lambda calcuation.
	//vectors to store skewness, standard deviation and mean of each gene/feature.
	vector< double > meanvec;
	
	//Objects to store the sums of variables that will be needed to perform linear regression.
	double thisSD, thisSkew;
	double SumSkew = 0;
	double SumSD = 0;
	double SumMean = 0;
	double SumSkewMean = 0;
	double SumSDMean = 0;
	double SumMeanSq = 0;
	
	double mean = 0;
	double M2 = 0;
	double M3 = 0;
	double delta, delta_n, term1;
	
	int_fast32_t i= 0, n = 0;
	double numData = 0;
	arma::mat::iterator it_end = GeneExp.end();
	//One pass linear regression with one pass variance, skewness
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		if (*it != 0) {
			delta = log1p(*it * lambda2) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			
			numData++;
		}
		n++;
		if (n == int(GeneExp.n_rows)) {
			//Here, calculate skewness and SD using 3rd and 2nd moments.
			thisSkew =  (sqrt(numData) * M3) / pow(M2,1.5);
			thisSD =  sqrt(M2/(numData - 1));
			
			//Linear regression
			meanvec.push_back(mean);
			SumSkew += thisSkew;
			SumSD += thisSD;
			SumMean += mean;
			SumSkewMean += thisSkew * mean;
			SumSDMean += thisSD * mean;
			SumMeanSq += pow(mean,2);
			
			mean = 0;
			M2 = 0;
			M3 = 0;
			numData = 0;
			n = 0;
			i++;
		}
	}
	//y = Mx + C, here are the M and C results of the linear regression
	double skewM = (SumSkewMean * GeneExp.n_cols - SumSkew * SumMean)/(SumMeanSq * GeneExp.n_cols - pow(SumMean,2) ) ;
	double skewC = (SumSkew - skewM * SumMean)/GeneExp.n_cols;
	double SDevM = (SumSDMean * GeneExp.n_cols - SumSD * SumMean)/(SumMeanSq * GeneExp.n_cols - pow(SumMean,2) );
	

	
	//integral of the linear equation of skewness
	double Skewzerointercept = -skewC/skewM;
	double Skewintegral;
	skewM = skewM/2;
	if (Skewzerointercept > meanvec.at(0) && Skewzerointercept < meanvec.at(GeneExp.n_cols-1)) {
		Skewintegral = (abs(skewM * (meanvec.at(GeneExp.n_cols-1) + Skewzerointercept) + skewC) * (meanvec.at(GeneExp.n_cols-1) - Skewzerointercept) + abs(skewM * (meanvec.at(0) + Skewzerointercept) + skewC)* (Skewzerointercept - meanvec.at(0)))/(meanvec.at(GeneExp.n_cols-1) - meanvec.at(0));
	} else {
		Skewintegral = abs(skewM * (meanvec.at(GeneExp.n_cols-1) + meanvec.at(0)) + skewC);
	}
	
	return pow(log1p(abs(SDevM))+1,2) + pow(log1p(Skewintegral)+1,2);
}
//This function outputs  F(lambda) and F(lambda + 1) together and in one pass, improving running speed.
static vector< double > SkewVar2 (arma::mat& GeneExp, const double& lambda2, double extra) {
	//double makes the program slower, but it is needed for accurate lambda calcuation.
	//vectors to store skewness, standard deviation and mean of each gene/feature.
	vector< double > meanvec;
	
	//Objects to store the sums of variables that will be needed to perform linear regression.
	double thisSD, thisSkew;
	double SumSkew = 0;
	double SumSD = 0;
	double SumMean = 0;
	double SumSkewMean = 0;
	double SumSDMean = 0;
	double SumMeanSq = 0;
	
	vector< double > meanvec_2;
	
	//Objects to store the sums of variables that will be needed to perform linear regression.
	double thisSD_2, thisSkew_2;
	double SumSkew_2 = 0;
	double SumSD_2 = 0;
	double SumMean_2 = 0;
	double SumSkewMean_2 = 0;
	double SumSDMean_2 = 0;
	double SumMeanSq_2 = 0;
	
	double mean = 0;
	double M2 = 0;
	double M3 = 0;
	double delta, delta_n, term1;
	
	double mean_2 = 0;
	double M2_2 = 0;
	double M3_2 = 0;
	double delta_2, delta_n_2, term1_2;
	int_fast32_t i= 0, n = 0;
	double numData = 0;
	arma::mat::iterator it_end = GeneExp.end();
	//One pass linear regression with one pass variance, skewness
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		if (*it != 0) {
			delta = log1p(*it * (lambda2 - extra + 1) ) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			
			delta_2 = log1p(*it * (lambda2 + extra) )- mean_2;
			delta_n_2 = delta_2 / (numData + 1);
			term1_2 = delta_2 * delta_n_2 * numData;
			mean_2 = mean_2 + delta_n_2;
			M3_2 = M3_2 + term1_2 * delta_n_2 * (numData - 1) - 3 * delta_n_2 * M2_2;
			M2_2 = M2_2 + term1_2;
			numData++;
		}
		n++;
		if (n == int(GeneExp.n_rows)) {
			//Here, calculate skewness and SD using 3rd and 2nd moments.
			thisSkew =  (sqrt(numData) * M3) / pow(M2,1.5);
			thisSD =  sqrt(M2/(numData - 1));
			
			//Linear regression
			meanvec.push_back(mean);
			SumSkew += thisSkew;
			SumSD += thisSD;
			SumMean += mean;
			SumSkewMean += thisSkew * mean;
			SumSDMean += thisSD * mean;
			SumMeanSq += pow(mean,2);
			
			//Here, calculate skewness and SD using 3rd and 2nd moments.
			thisSkew_2 =  (sqrt(numData) * M3_2) / pow(M2_2,1.5);
			thisSD_2 =  sqrt(M2_2/(numData - 1));
			
			//Linear regression
			meanvec_2.push_back( mean_2 );
			SumSkew_2 += thisSkew_2;
			SumSD_2 += thisSD_2;
			SumMean_2 += mean_2;
			SumSkewMean_2 += thisSkew_2 * mean_2;
			SumSDMean_2 += thisSD_2 * mean_2;
			SumMeanSq_2 += pow(mean_2,2);
			mean = 0;
			M2 = 0;
			M3 = 0;

			mean_2 = 0;
			M2_2 = 0;
			M3_2 = 0;
			
			numData = 0;
			n = 0;
			i++;
		}
	}
	//y = Mx + C, here are the M and C results of the linear regression
	double skewM = (SumSkewMean * GeneExp.n_cols - SumSkew * SumMean)/(SumMeanSq * GeneExp.n_cols - pow(SumMean,2) ) ;
	double skewC = (SumSkew - skewM * SumMean)/GeneExp.n_cols;
	double SDevM = (SumSDMean * GeneExp.n_cols - SumSD * SumMean)/(SumMeanSq * GeneExp.n_cols - pow(SumMean,2) );
	
	double skewM_2 = (SumSkewMean_2 * GeneExp.n_cols - SumSkew_2 * SumMean_2)/(SumMeanSq_2 * GeneExp.n_cols - pow(SumMean_2,2) ) ;
	double skewC_2 = (SumSkew_2 - skewM_2 * SumMean_2)/(GeneExp.n_cols);
	double SDevM_2 = (SumSDMean_2 * GeneExp.n_cols - SumSD_2 * SumMean_2)/(SumMeanSq_2 * GeneExp.n_cols - pow(SumMean_2,2) );
	
	
	//integral of the linear equation of skewness
	double Skewzerointercept = -skewC/skewM;
	double Skewintegral;
	skewM = skewM/2;
	if (Skewzerointercept > meanvec.at(0) && Skewzerointercept < meanvec.at(GeneExp.n_cols-1)) {
		Skewintegral = (abs(skewM * (meanvec.at(GeneExp.n_cols-1) + Skewzerointercept) + skewC) * (meanvec.at(GeneExp.n_cols-1) - Skewzerointercept) + abs(skewM * (meanvec.at(0) + Skewzerointercept) + skewC)* (Skewzerointercept - meanvec.at(0)))/(meanvec.at(GeneExp.n_cols-1) - meanvec.at(0));
	} else {
		Skewintegral = abs(skewM * (meanvec.at(GeneExp.n_cols-1) + meanvec.at(0)) + skewC);
	}
	
	//integral of the linear equation of skewness
	double Skewzerointercept_2 = -skewC_2/skewM_2;
	double Skewintegral_2;
	skewM_2 = skewM_2/2;
	if (Skewzerointercept_2 > meanvec_2.at(0) && Skewzerointercept_2 < meanvec_2.at(GeneExp.n_cols-1)) {
		Skewintegral_2 = (abs(skewM_2 * (meanvec_2.at(GeneExp.n_cols-1) + Skewzerointercept_2) + skewC_2) * (meanvec_2.at(GeneExp.n_cols-1) - Skewzerointercept_2) + abs(skewM_2 * (meanvec_2.at(0) + Skewzerointercept_2) + skewC_2)* (Skewzerointercept_2 - meanvec_2.at(0)))/(meanvec_2.at(GeneExp.n_cols-1) - meanvec_2.at(0));
	} else {
		Skewintegral_2 = abs(skewM_2 * (meanvec_2.at(GeneExp.n_cols-1) + meanvec_2.at(0)) + skewC_2);
	}

	vector< double > Answer;
	Answer.push_back( pow(log1p(abs(SDevM))+1,2) + pow(log1p(Skewintegral)+1,2) );
	Answer.push_back( pow(log1p(abs(SDevM_2))+1,2) + pow(log1p(Skewintegral_2)+1,2) );
	return Answer;
}
SEXP SkewVar2Cpp(SEXP xSEXP,SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double lambda2 = Rcpp::as<double>(ySEXP);
	return Rcpp::wrap(SkewVar2(GeneExp, lambda2, 1));
}

//Legacy: This function output S(lambda) and V(lambda) functions for testing. 
static arma::vec SkewAVar (const arma::mat& GeneExp, const double& lambda2) {
	//vectors to store skewness, standard deviation and mean of each gene/feature.
	double lambda = (double) lambda2;
	
	double thisSD, thisSkew;
	vector< double > meanvec;
	//Objects to store the sums of variables that will be needed to perform linear regression.
	double SumSkew = 0;
	double SumSD = 0;
	double SumMean = 0;
	double SumSkewMean = 0;
	double SumSDMean = 0;
	double SumMeanSq = 0;
	
	double mean = 0;
	double M2 = 0;
	double M3 = 0;
	double numData = 0;
	double delta, delta_n, term1;
		
	//One pass linear regression with one pass variance, skewness
	for (int_fast32_t i = 0; i < int_fast32_t(GeneExp.n_rows); i ++) {
		mean = 0;
		M2 = 0;
		M3 = 0;
		numData = 0;
		for (int_fast32_t n = 0; n < int_fast32_t(GeneExp.n_cols); n++) {
			if (GeneExp.at(i,n) != 0) {
				delta = log1p(GeneExp.at(i,n) * lambda) - mean;
				delta_n = delta / (numData + 1);
				term1 = delta * delta_n * numData;
				mean = mean + delta_n;
				M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
				M2 = M2 + term1;
				numData++;
			}
		}
		//Here, calculate skewness and SD using 3rd and 2nd moments.
		thisSkew =  (sqrt(numData) * M3) / pow(M2,1.5);
		thisSD = sqrt(M2/(numData - 1));
		
		//Linear regression
		meanvec.push_back(mean);
		SumSkew += thisSkew;
		SumSD += thisSD;
		SumMean += meanvec.at(i);
		SumSkewMean += thisSkew * mean;
		SumSDMean += thisSD * mean;
		SumMeanSq += pow(meanvec.at(i),2);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	double skewM = (SumSkewMean * GeneExp.n_rows - SumSkew * SumMean)/(SumMeanSq * GeneExp.n_rows - pow(SumMean,2) ) ;
	double skewC = (SumSkew - skewM * SumMean)/GeneExp.n_rows;
	double SDevM = (SumSDMean * GeneExp.n_rows - SumSD * SumMean)/(SumMeanSq * GeneExp.n_rows - pow(SumMean,2) );
	
	//integral of the linear equation of skewness
	double Skewzerointercept = -skewC/skewM;
	double Skewintegral;
	skewM = skewM/2;
	if (Skewzerointercept > meanvec.at(0) && Skewzerointercept < meanvec.at(GeneExp.n_rows-1)) {
		Skewintegral = (abs(skewM * (meanvec.at(GeneExp.n_rows-1) + Skewzerointercept) + skewC) * (meanvec.at(GeneExp.n_rows-1) - Skewzerointercept) + abs(skewM * (meanvec.at(0) + Skewzerointercept) + skewC)* (Skewzerointercept - meanvec.at(0)))/(meanvec.at(GeneExp.n_rows-1) - meanvec.at(0));
	} else {
		Skewintegral = abs(skewM * (meanvec.at(GeneExp.n_rows-1) + meanvec.at(0)) + skewC) ;
	}
	arma::vec Answer(2);
	Answer.at(0) = SDevM;
	Answer.at(1) = Skewintegral;
	return Answer;
}
SEXP SkewAVarCpp(SEXP xSEXP,SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double lambda2 = Rcpp::as<double>(ySEXP);
	return Rcpp::wrap(SkewAVar(GeneExp, lambda2));
}

//Given a range of lambda, this function finds lambda that minimizes F(lambda) (see article) based on the expression matrix by using a history recording binary search.
static double LocalSearch(arma::mat& GeneExp, double minBound, double maxBound, double& smallest, double search_exponent) {
	minBound = round(minBound);
	maxBound = round(maxBound);
	//cout << "here1 " << minBound<<   " " << maxBound << endl;
	double GminBound = minBound, GmaxBound = maxBound;
	double midBound = round((minBound + maxBound)/2);
	
	//Initialize smallestBound. This object is for the program to remember the smallest F(lambda) ever calculated. If the local minimal found is larger that this smallestBound, the boundaries will be reset and binary search will be rerun using a smaller boundary, with the smallestBound as center.
	double smallestBound;
	double extra = 1;
	vector< double > om;
	while (true) {
		om = SkewVar2(GeneExp, midBound, extra);
		if (om.at(0) > om.at(1)) {
			minBound = midBound;
			smallest = om.at(1);
			smallestBound = round(midBound + extra);
			break;
		} else if (om.at(0) < om.at(1)) {
			maxBound = midBound;
			smallest = om.at(0);
			smallestBound = round(midBound - extra + 1);
			break;
		}
		extra++;
	}
	midBound = round((minBound + maxBound)/2);
	int_fast32_t index = 0;
	while (smallestBound != midBound) {
		//cout << "here1 " << smallestBound<<   " " << midBound << endl;
		if (index > 0) {
			//Reset boundary to center at smallestBound
			minBound = round(smallestBound - (smallestBound - GminBound)/pow(search_exponent,index));
			maxBound = round(smallestBound + (GmaxBound - smallestBound)/pow(search_exponent,index));
			midBound = round((minBound + maxBound)/2);
			//cout << minBound << " " << maxBound << endl;
		}
		index++;
		//Binary search for smallest F(lambda)
		while (maxBound > midBound) {
			extra = 1;
			while (true) {
				om = SkewVar2(GeneExp, midBound, extra);
				double thisom = om.at(0);
				if (om.at(0) > om.at(1)) {
					minBound = midBound;
					if (smallest > om.at(1)) {
						smallest = om.at(1);
						smallestBound = round(midBound + extra);
					}
					break;
				} else if (om.at(0) < om.at(1)) {
					maxBound = midBound;
					if (smallest > om.at(0)) {
						smallest = om.at(0);
						smallestBound = round(midBound - extra + 1);
					}
					break;
				}
				if (thisom == smallest) {
					smallestBound = midBound;
				}
				extra++;
			}
			midBound = round((minBound + maxBound)/2);
		}
	}
	return midBound;
}
static double LocalSearch_legacy(const arma::mat& GeneExp, double minBound, double maxBound, double& smallest, double search_exponent) {
	minBound = round(minBound);
	maxBound = round(maxBound);
	//cout << "here1 " << minBound<<   " " << maxBound << endl;
	double GminBound = minBound, GmaxBound = maxBound;
	double midBound = round((minBound + maxBound)/2);
	double om = SkewVar(GeneExp, midBound);
	double om2 = SkewVar(GeneExp,(midBound + 1));
	//Initialize smallestBound. This object is for the program to remember the smallest F(lambda) ever calculated. If the local minimal found is larger that this smallestBound, the boundaries will be reset and binary search will be rerun using a smaller boundary, with the smallestBound as center.
	double smallestBound;
	if (om > om2) {
		minBound = midBound + 1;
		smallest = om2;
		smallestBound = midBound + 1;
	} else {
		maxBound = midBound;
		smallest = om;
		smallestBound = midBound;
	}
	midBound = round((minBound + maxBound)/2);
	int_fast32_t index = 0;
	while (smallestBound != midBound) {
		//cout << "here1 " << smallestBound<<   " " << midBound << endl;
		if (index > 0) {
			//Reset boundary to center at smallestBound
			minBound = round(smallestBound - (smallestBound - GminBound)/pow(search_exponent,index));
			maxBound = round(smallestBound + (GmaxBound - smallestBound)/pow(search_exponent,index));
			midBound = round((minBound + maxBound)/2);
			//cout << minBound << " " << maxBound << endl;
		}
		index++;
		//Binary search for smallest F(lambda)
		while (maxBound > midBound) {
			om = SkewVar(GeneExp, midBound);
			om2 = SkewVar(GeneExp,midBound + 1);
			if (om > om2) {
				minBound = midBound + 1;
				if (om2 < smallest) {
					smallest = om2;
					smallestBound = round(midBound + 1);
				}
			} else {
				maxBound = midBound;
				if (om < smallest) {
					smallest = om;
					smallestBound = midBound;
				}
			}
			midBound = round((minBound + maxBound)/2);
		}
	}
	return smallestBound;
}

//Using the LocalSearch function, performs iterated local search to find minimal lambda.
SEXP LocateLambdaCpp(SEXP xSEXP,SEXP ySEXP, SEXP zSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double search_exponent = Rcpp::as<double>(ySEXP);
	double maxBound = Rcpp::as<double>(zSEXP);
   //Boundary of ILS are minBound and maxBound
	//Here, we define the range of lambda based on the dataset.
	double OriginalmaxBound = maxBound;
	double minBound = 1;
	//Find local minima
	double localminIntegral;
	//Intital minimal
	double localmin = LocalSearch(GeneExp, minBound, maxBound,localminIntegral,search_exponent);
	
	//First, iterated local search starting on the left hand side of the local minima.  
	//LHS search
	int_fast32_t searchIndex = 2;
	//1.perturbing the current local minimum;
	double newminBound = localmin -  10 * pow(search_exponent,searchIndex);
	if (newminBound < 1) {
		newminBound = minBound;
	}
	double lastnewminBound = localmin;
	double newmin;
	double newminIntegral;
	double finallocalmin = localmin + 1;
	int_fast32_t numMaxBoundIncrease = 0;
	while ((finallocalmin != localmin || maxBound - 1 <= localmin || minBound + 1 >= localmin) && numMaxBoundIncrease <= 3) {
		if ((minBound + 1 >= localmin || maxBound - 1 <= localmin) && newminBound < minBound) {
			//enlarge maxBound if global minimum is found near the boundaries.
			maxBound = maxBound * 10;
			localmin = round((minBound + maxBound)/2);
			lastnewminBound = localmin;
			localminIntegral = SkewVar1(GeneExp, localmin);
			searchIndex = 2;
			newminBound = localmin -  10 * pow(search_exponent,searchIndex);
			if (newminBound < minBound) {
				newminBound = minBound;
			}
			numMaxBoundIncrease++;
			//cout << "maxBound enlarged to " << maxBound << endl;
		} else if (newminBound < minBound) {
			searchIndex = 2;
			newminBound = localmin -  10 * pow(search_exponent,searchIndex);
			lastnewminBound = localmin;
			if (newminBound < minBound) {
				newminBound = minBound;
			}
		}
		finallocalmin = localmin;
		//Note: instead of randomly choosing a new range to search, we perturb from the local minimal by increasing the range exponentially. We first search through left hand side(LHS), then search through right hand side(RHS) of the local minimal. The stop condition is that no new local minimal is found after a full LHS and RHS search. This allow us to keep the running time around mlog(n). We use a history recording binary search algorithm to find local minimum within a range.
		while (newminBound >= minBound) {
			//cout << "LHS " << searchIndex <<  " " << newmin <<endl;
			//2.applying local search after starting from the modified solution.
			newmin = LocalSearch(GeneExp, newminBound, lastnewminBound, newminIntegral,search_exponent);
			//If new minimum is smaller than the local minimal, reset local minimal and searchIndex.
			if (newminIntegral < localminIntegral) {
				localmin = newmin;
				localminIntegral = newminIntegral;
				lastnewminBound = (newminBound + lastnewminBound)/2;
				searchIndex++;
				//1.perturbing the current local minimum;
				if (newminBound != minBound) {
					newminBound = localmin -  10 * pow(search_exponent,searchIndex);
					if (newminBound < minBound) {
						newminBound = minBound;
					}
				} else {
					newminBound = minBound - 1;
				}
			} else {
				searchIndex++;
				lastnewminBound = (newminBound + lastnewminBound)/2;
				//1.perturbing the current local minimum;
				if (newminBound != minBound) {
					newminBound = localmin -  10 * pow(search_exponent,searchIndex);
					if (newminBound < minBound) {
						newminBound = minBound;
					}
				} else {
					newminBound = minBound - 1;
				}
			}
		}

		
		//Lastly, iterated local search on the right hand side of the local minima.  
		//RHS search
		searchIndex = 2;
		//1.perturbing the current local maximum;
		double newmaxBound = localmin + 10 * pow(search_exponent,searchIndex);
		if (newmaxBound > maxBound) {
			newmaxBound = maxBound;
		}
		double lastnewmaxBound = localmin;
		double newmax;
		double newmaxIntegral;
		while (newmaxBound <= maxBound) {
			//cout << "RHS " << searchIndex <<  " " << newmax <<endl;
			//2.applying local search after starting from the modified solution.
			newmax = LocalSearch(GeneExp, lastnewmaxBound, newmaxBound, newmaxIntegral,search_exponent);
			//If new maximum is smaller than the local maximal, reset local maximal and searchIndex.
			if (newmaxIntegral < localminIntegral) {
				localmin = newmax;
				localminIntegral = newmaxIntegral;
				lastnewmaxBound = (newmaxBound + lastnewmaxBound)/2;
				searchIndex++;
				//1.perturbing the current local maximum;
				if (newmaxBound != maxBound) {
					newmaxBound = localmin +  10 * pow(search_exponent,searchIndex);
					if (newmaxBound > maxBound) {
						
						newmaxBound = maxBound;
					}
				} else {
					newmaxBound = maxBound + 1;
				}
			} else {
				searchIndex++;
				lastnewmaxBound = (newmaxBound + lastnewmaxBound)/2;
				//1.perturbing the current local maximum;
				if (newmaxBound != maxBound) {
					newmaxBound = localmin +  10 * pow(search_exponent,searchIndex);
					if (newmaxBound > maxBound) {
						newmaxBound = maxBound;
					}
				} else {
					newmaxBound = maxBound + 1;
				}
			}
		}
	}
	if (numMaxBoundIncrease == 4 && (maxBound - 1 <= localmin || minBound + 1 >= localmin)) {
		localmin = OriginalmaxBound;
	}
	//cout << "LMI " << localminIntegral << endl;
	return Rcpp::wrap(localmin);
}
SEXP LocateLambdaCpp_legacy(SEXP xSEXP,SEXP ySEXP, SEXP zSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double search_exponent = Rcpp::as<double>(ySEXP);
	double maxBound = Rcpp::as<double>(zSEXP);
   //Boundary of ILS are minBound and maxBound
	//Here, we define the range of lambda based on the dataset.
	double OriginalmaxBound = maxBound;
	double minBound = 1;
	//Find local minima
	double localminIntegral;
	//Intital minimal
	double localmin = LocalSearch_legacy(GeneExp, minBound, maxBound,localminIntegral,search_exponent);
	
	//First, iterated local search starting on the left hand side of the local minima.  
	//LHS search
	int_fast32_t searchIndex = 2;
	//1.perturbing the current local minimum;
	double newminBound = localmin -  10 * pow(search_exponent,searchIndex);
	if (newminBound < 1) {
		newminBound = minBound;
	}
	double lastnewminBound = localmin;
	double newmin;
	double newminIntegral;
	double finallocalmin = localmin + 1;
	int_fast32_t numMaxBoundIncrease = 0;
	while ((finallocalmin != localmin || maxBound - 1 <= localmin || minBound + 1 >= localmin) && numMaxBoundIncrease <= 5) {
		if ((minBound + 1 >= localmin || maxBound - 1 <= localmin) && newminBound < minBound) {
			//enlarge maxBound if global minimum is found near the boundaries.
			maxBound = maxBound * 10;
			localmin = round((minBound + maxBound)/2);
			lastnewminBound = localmin;
			localminIntegral = SkewVar(GeneExp, localmin);
			searchIndex = 2;
			newminBound = localmin -  10 * pow(search_exponent,searchIndex);
			if (newminBound < minBound) {
				newminBound = minBound;
			}
			numMaxBoundIncrease++;
			//cout << "maxBound enlarged to " << maxBound << endl;
		} else if (newminBound < minBound) {
			searchIndex = 2;
			newminBound = localmin -  10 * pow(search_exponent,searchIndex);
			lastnewminBound = localmin;
			if (newminBound < minBound) {
				newminBound = minBound;
			}
		}
		finallocalmin = localmin;
		//Note: instead of randomly choosing a new range to search, we perturb from the local minimal by increasing the range exponentially. We first search through left hand side(LHS), then search through right hand side(RHS) of the local minimal. The stop condition is that no new local minimal is found after a full LHS and RHS search. This allow us to keep the running time around mlog(n). We use a history recording binary search algorithm to find local minimum within a range.
		while (newminBound >= minBound) {
			//cout << "LHS " << searchIndex <<  " " << newmin <<endl;
			//2.applying local search after starting from the modified solution.
			newmin = LocalSearch_legacy(GeneExp, newminBound, lastnewminBound, newminIntegral,search_exponent);
			//If new minimum is smaller than the local minimal, reset local minimal and searchIndex.
			if (newminIntegral < localminIntegral) {
				localmin = newmin;
				localminIntegral = newminIntegral;
				lastnewminBound = (newminBound + lastnewminBound)/2;
				searchIndex++;
				//1.perturbing the current local minimum;
				if (newminBound != minBound) {
					newminBound = localmin -  10 * pow(search_exponent,searchIndex);
					if (newminBound < minBound) {
						newminBound = minBound;
					}
				} else {
					newminBound = minBound - 1;
				}
			} else {
				searchIndex++;
				lastnewminBound = (newminBound + lastnewminBound)/2;
				//1.perturbing the current local minimum;
				if (newminBound != minBound) {
					newminBound = localmin -  10 * pow(search_exponent,searchIndex);
					if (newminBound < minBound) {
						newminBound = minBound;
					}
				} else {
					newminBound = minBound - 1;
				}
			}
		}

		
		//Lastly, iterated local search on the right hand side of the local minima.  
		//RHS search
		searchIndex = 2;
		//1.perturbing the current local maximum;
		double newmaxBound = localmin + 10 * pow(search_exponent,searchIndex);
		if (newmaxBound > maxBound) {
			newmaxBound = maxBound;
		}
		double lastnewmaxBound = localmin;
		double newmax;
		double newmaxIntegral;
		while (newmaxBound <= maxBound) {
			//cout << "RHS " << searchIndex <<  " " << newmax <<endl;
			//2.applying local search after starting from the modified solution.
			newmax = LocalSearch_legacy(GeneExp, lastnewmaxBound, newmaxBound, newmaxIntegral,search_exponent);
			//If new maximum is smaller than the local maximal, reset local maximal and searchIndex.
			if (newmaxIntegral < localminIntegral) {
				localmin = newmax;
				localminIntegral = newmaxIntegral;
				lastnewmaxBound = (newmaxBound + lastnewmaxBound)/2;
				searchIndex++;
				//1.perturbing the current local maximum;
				if (newmaxBound != maxBound) {
					newmaxBound = localmin +  10 * pow(search_exponent,searchIndex);
					if (newmaxBound > maxBound) {
						
						newmaxBound = maxBound;
					}
				} else {
					newmaxBound = maxBound + 1;
				}
			} else {
				searchIndex++;
				lastnewmaxBound = (newmaxBound + lastnewmaxBound)/2;
				//1.perturbing the current local maximum;
				if (newmaxBound != maxBound) {
					newmaxBound = localmin +  10 * pow(search_exponent,searchIndex);
					if (newmaxBound > maxBound) {
						newmaxBound = maxBound;
					}
				} else {
					newmaxBound = maxBound + 1;
				}
			}
		}
	}
	if (numMaxBoundIncrease == 6 && (maxBound - 1 <= localmin || minBound + 1 >= localmin)) {
		localmin = OriginalmaxBound;
	}
	//cout << "LMI " << localminIntegral << endl;
	return Rcpp::wrap(localmin);
}

//Get Slope from x and y vectors
SEXP getSlopeCpp(SEXP xSEXP, SEXP ySEXP) {
	arma::vec xvec = Rcpp::as<arma::vec>(xSEXP);
	arma::vec yvec = Rcpp::as<arma::vec>(ySEXP);
	
	double SumX = 0;
	double SumXSq = 0;
	double SumXY = 0;
	double SumY = 0;
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int_fast32_t i = 0; i < int_fast32_t(xvec.n_elem); i ++) {
		SumX += xvec.at(i);
		SumXSq += pow(xvec.at(i),2);
		SumXY += xvec.at(i) * yvec.at(i);
		SumY += yvec.at(i);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	double Slope = (xvec.n_elem * SumXY - SumY * SumX)/(xvec.n_elem * SumXSq - pow(SumX,2) ) ;
	//double constant = (SumY.at(n) - Slope * SumX)/GeneExp.n_rows;
	return Rcpp::wrap(Slope);
}

SEXP LinearRegressionCpp(SEXP xSEXP, SEXP ySEXP) {
	arma::vec xvec = Rcpp::as<arma::vec>(xSEXP);
	arma::vec yvec = Rcpp::as<arma::vec>(ySEXP);
	
	double SumX = 0;
	double SumXSq = 0;
	double SumXY = 0;
	double SumY = 0;
	
	//One pass linear regression
	for (int_fast32_t i = 0; i < int_fast32_t(xvec.n_elem); i++) {
		SumX += xvec.at(i);
		SumXSq += pow(xvec.at(i),2);
		SumXY += xvec.at(i) * yvec.at(i);
		SumY += yvec.at(i);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	arma::vec Answer(2);
	Answer.at(0) = (xvec.n_elem * SumXY - SumY * SumX)/(xvec.n_elem * SumXSq - pow(SumX,2) ) ;
	Answer.at(1) = (SumY - Answer.at(0) * SumX)/xvec.n_elem;
	return Rcpp::wrap(Answer);
}
SEXP LinearRegressionZeroCpp(SEXP xSEXP, SEXP ySEXP) {
	arma::vec xvec = Rcpp::as<arma::vec>(xSEXP);
	arma::vec yvec = Rcpp::as<arma::vec>(ySEXP);
	
	double SumXSq = 0;
	double SumXY = 0;
	
	//One pass linear regression
	for (int_fast32_t i = 0; i < int_fast32_t(xvec.n_elem); i ++) {
		SumXY += xvec.at(i) *  yvec.at(i);
		SumXSq += pow(xvec.at(i),2);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	arma::vec Answer(2);
	Answer.at(0) = 0;
	Answer.at(1) = SumXY/SumXSq;
	return Rcpp::wrap(Answer);
}
SEXP LinearRegressionFPCpp(SEXP xSEXP, SEXP ySEXP,SEXP x1SEXP, SEXP y1SEXP) {
	arma::vec xvec = Rcpp::as<arma::vec>(xSEXP);
	arma::vec yvec = Rcpp::as<arma::vec>(ySEXP);
	
	double x1 = Rcpp::as<double>(x1SEXP);
	double y1 = Rcpp::as<double>(y1SEXP);
	
	double SumX = 0;
	double SumXSq = 0;
	double SumXY = 0;
	double SumY = 0;

	//One pass linear regression
	for (int_fast32_t i = 0; i < int_fast32_t(xvec.n_elem); i ++) {
		SumX += xvec.at(i);
		SumXSq += pow(xvec.at(i),2);
		SumXY += xvec.at(i) * yvec.at(i);
		SumY += yvec.at(i);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	arma::vec Answer(2);
	Answer.at(1) = (x1*y1*SumX-y1*SumXSq-pow(x1,2)*SumY+x1*SumXY)/(2*x1*SumX-SumXSq-pow(x1,2)*xvec.n_elem);
	Answer.at(0) = (y1 - Answer.at(1))/x1;
	return Rcpp::wrap(Answer);
}

//Algorithm to normalize batch effect
SEXP BatchEffectCpp(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP z2SEXP){
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	arma::mat Output = Rcpp::as<arma::mat>(ySEXP);
	arma::vec meanvec = Rcpp::as<arma::vec>(zSEXP);
	double BE_strength = Rcpp::as<double>(z2SEXP);
	
	//Objects to store the sums of variables that will be needed to perform linear regression.

	arma::vec SumX(GeneExp.n_rows);
	SumX.fill(0);
	arma::vec SumXSq(GeneExp.n_rows);
	SumXSq.fill(0);
	arma::vec SumXY(GeneExp.n_rows);
	SumXY.fill(0);
	arma::vec SumY(GeneExp.n_rows);
	SumY.fill(0);
	arma::vec nRows(GeneExp.n_rows);
	nRows.fill(0);
	
	double store;
	arma::mat::iterator it_end = GeneExp.end();
	int_fast32_t i=0, n=0;
	//One pass linear regression on all rows
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		if (*it != 0) {
			store = log(*it);
			SumXY.at(n) += store *  meanvec.at(i);
			SumX.at(n) += store;
			nRows.at(n)++;
			SumXSq.at(n) += pow( store,2);
			SumY.at(n) += meanvec.at(i);
		}
		n++;
		if (n == int(GeneExp.n_rows)) {
			n=0;
			i++;
		}
	}
	
	//y = Mx + C, here are the M and C results of the linear regression
	double add = (1/BE_strength) - 1;
	arma::vec Slope(GeneExp.n_rows);
	arma::vec constant(GeneExp.n_rows);
	for (int_fast32_t k = 0; k < int_fast32_t(GeneExp.n_rows); k++) {
		Slope.at(k) = ((nRows.at(k) * SumXY.at(k) - SumY.at(k) * SumX.at(k))/(nRows.at(k) * SumXSq.at(k) - pow(SumX.at(k),2) )+ add) * BE_strength;
		constant.at(k) = ((SumY.at(k) - Slope.at(k) * SumX.at(k))/nRows.at(k)) * BE_strength;
	}
	
	
	arma::mat::iterator it_end2 = Output.end();
	i=0;
	n=0;
	//Update Output matrix for output
	for (arma::mat::iterator it2 = Output.begin(); it2 != it_end2; ++it2) {
		//std::cout << (*it) << std::endl;
		if (*it2 != 0) {
			*it2 = exp(Slope.at(n) * log(*it2) + constant.at(n));
		}
		n++;
		if (n == int(Output.n_rows)) {
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(Output);
}

//Misc R functions for mean, sd and skewness calculations. Different ones are used in different situations. They are all implemented to improve running time.
SEXP colVarsCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::vec Answer(GeneExp.n_cols);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	int_fast32_t i=0, n=0;
	arma::mat::iterator it_end = GeneExp.end();
	double mean = 0;
	double M2 = 0;
	double delta, delta_n, term1;
	int_fast32_t numData = 0;
	//One pass
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		if (*it != 0) {
			numData++;
		}
		delta = *it - mean;
		delta_n = delta / (n + 1);
		term1 = delta * delta_n * n;
		mean = mean + delta_n;
		M2 = M2 + term1;
		n++;
		if (n == int(GeneExp.n_rows)) {
			switch (numData) {
				case 0:
					Answer.at(i) = std::numeric_limits<double>::quiet_NaN();
					break;
				case 1:
					Answer.at(i) = std::numeric_limits<double>::quiet_NaN();
					break;
				default:
					Answer.at(i) = M2/(GeneExp.n_rows - 1);
			}
			numData = 0;
			M2 = 0;
			mean = 0;
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP colSDsCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::vec Answer(GeneExp.n_cols);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	int_fast32_t i=0, n=0;
	arma::mat::iterator it_end = GeneExp.end();
	double mean = 0;
	double M2 = 0;
	double delta, delta_n, term1;
	int_fast32_t numData = 0;
	//One pass
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		if (*it != 0) {
			numData++;
		}
		delta = *it - mean;
		delta_n = delta / (n + 1);
		term1 = delta * delta_n * n;
		mean = mean + delta_n;
		M2 = M2 + term1;
		n++;
		if (n == int(GeneExp.n_rows)) {
			switch (numData) {
				case 0:
					Answer.at(i) = std::numeric_limits<double>::quiet_NaN();
					break;
				case 1:
					Answer.at(i) = std::numeric_limits<double>::quiet_NaN();
					break;
				default:
					Answer.at(i) = sqrt(M2/(GeneExp.n_rows - 1));
			}
			numData = 0;
			M2 = 0;
			mean = 0;
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP colMeanSDCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(2,GeneExp.n_cols);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	int_fast32_t i=0, n=0;
	arma::mat::iterator it_end = GeneExp.end();
	double mean = 0;
	double M2 = 0;
	double delta, delta_n, term1;
	int_fast32_t numData = 0;
	//One pass
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		if (*it != 0) {
			numData++;
		}
		delta = *it - mean;
		delta_n = delta / (n + 1);
		term1 = delta * delta_n * n;
		mean = mean + delta_n;
		M2 = M2 + term1;
		n++;
		if (n == int(GeneExp.n_rows)) {
			switch (numData) {
				case 0:
					Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
					Answer.at(0,i) = 0;
					break;
				case 1:
					Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
					Answer.at(0,i) = mean;
					break;
				default:
					Answer.at(1,i) = sqrt(M2/(GeneExp.n_rows - 1));
					Answer.at(0,i) = mean;
			}
			numData = 0;
			M2 = 0;
			mean = 0;
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(Answer);
}

SEXP NZcolLog1pMeanSDCpp (SEXP xSEXP, SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double lambda = Rcpp::as<double>(ySEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(2,GeneExp.n_cols);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	int_fast32_t i=0, n=0;
	arma::mat::iterator it_end = GeneExp.end();
	double mean = 0;
	double M2 = 0;
	double delta, delta_n, term1;
	int numdata = 0;
	//One pass
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		if (*it != 0) {
			delta = log1p(*it * lambda) - mean;
			delta_n = delta / (n + 1);
			term1 = delta * delta_n * n;
			mean = mean + delta_n;
			M2 = M2 + term1;
			numdata++;
		}
		n++;
		if (n == int(GeneExp.n_rows)) {
			Answer.at(1,i) = sqrt(M2/(numdata - 1));
			Answer.at(0,i) = mean;
			numdata = 0;
			M2 = 0;
			mean = 0;
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(Answer);
}

SEXP NZcolLogMeanSDSkewCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(3,GeneExp.n_cols);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	int_fast32_t i=0, n=0;
	arma::mat::iterator it_end = GeneExp.end();
	double mean = 0;
	double M2 = 0;
	double M3 = 0;
	double delta, delta_n, term1;
	int_fast32_t numData = 0;
	//One pass
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		if (*it != 0) {
			delta = log(*it) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		n++;
		if (n == int(GeneExp.n_rows)) {
			switch (numData) {
				case 0:
					Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
					Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
					Answer.at(0,i) = 0;
					break;
				case 1:
					Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
					Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
					Answer.at(0,i) = mean;
					break;
				case 2:
					Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
					Answer.at(1,i) = sqrt(M2/(GeneExp.n_rows - 1));
					Answer.at(0,i) = mean;
					break;
				default:
					Answer.at(2,i) = (sqrt(GeneExp.n_rows) * M3) / pow(M2,1.5);
					Answer.at(1,i) = sqrt(M2/(GeneExp.n_rows - 1));
					Answer.at(0,i) = mean;
			}
			numData = 0;
			mean = 0;
			M2 = 0;
			M3 = 0;
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP NZcolMeansCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::vec Answer(GeneExp.n_cols);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	int_fast32_t i=0, n=0;
	arma::mat::iterator it_end = GeneExp.end();
	double mean = 0;
	double delta, delta_n;
	int_fast32_t numData = 0;
	//One pass
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		if (*it != 0) {
			delta = *it - mean;
			delta_n = delta / (numData + 1);
			mean = mean + delta_n;
			numData++;
		}
		n++;
		if (n == int(GeneExp.n_rows)) {
			Answer.at(i) = mean;
			numData = 0;
			mean = 0;
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP NZcolMeanSDCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(2,GeneExp.n_cols);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	int_fast32_t i=0, n=0;
	arma::mat::iterator it_end = GeneExp.end();
	double mean = 0;
	double M2 = 0;
	double delta, delta_n, term1;
	int_fast32_t numData = 0;
	//One pass
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		if (*it != 0) {
			delta = *it - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M2 = M2 + term1;
			numData++;
		}
		n++;
		if (n == int(GeneExp.n_rows)) {
			switch (numData) {
				case 0:
					Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
					Answer.at(0,i) = 0;
					break;
				case 1:
					Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
					Answer.at(0,i) = mean;
					break;
				default:
					Answer.at(1,i) = sqrt(M2/(GeneExp.n_rows - 1));
					Answer.at(0,i) = mean;
			}
			numData = 0;
			M2 = 0;
			mean = 0;
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(Answer);
}

SEXP WNZcolMeansCpp(SEXP xSEXP, SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	
	arma::vec weights = Rcpp::as<arma::vec>(ySEXP);
	arma::vec Answer(GeneExp.n_cols);
	
	int_fast32_t i=0, n=0;
	arma::mat::iterator it_end = GeneExp.end();
	double mean = 0, sumweight = 0;
	double delta, delta_n;
	int_fast32_t numData = 0;
	//One pass
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		if (*it != 0) {
			delta = *it * weights.at(n) - mean;
			delta_n = delta / (numData + 1);
			mean = mean + delta_n;
			numData++;
			sumweight += weights.at(n);
		}
		n++;
		if (n == int(GeneExp.n_rows)) {
			if (numData > 0) {
				Answer.at(i) = (mean * numData)/(sumweight);
			} else {
				Answer.at(i) = 0;
			}
			
			
			numData = 0;
			mean = 0;
			sumweight = 0;
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP WcolMeansCpp(SEXP xSEXP, SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	
	arma::vec weights = Rcpp::as<arma::vec>(ySEXP);
	double sumweight = accu(weights);
	arma::vec Answer(GeneExp.n_cols);
	
	int_fast32_t i=0, n=0;
	arma::mat::iterator it_end = GeneExp.end();
	double mean = 0;
	double delta, delta_n;
	//One pass
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		delta = *it * weights.at(n) - mean;
		delta_n = delta / (n + 1);
		mean = mean + delta_n;
		n++;
		if (n == int(GeneExp.n_rows)) {
			Answer.at(i) = (mean * n)/(sumweight);
			
			mean = 0;
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(Answer);
}

SEXP popZScoreCPP(SEXP xSEXP) {
	arma::vec DataVector = Rcpp::as<arma::vec>(xSEXP);
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	double mean = 0;
	double M2 = 0;
	double delta, delta_n, term1;
	for (int_fast32_t n = 0; n < int_fast32_t(DataVector.n_elem); n++) {
		delta = DataVector.at(n) - mean;
		delta_n = delta / (n + 1);
		term1 = delta * delta_n * n;
		mean = mean + delta_n;
		M2 = M2 + term1;
	}
	double popSD =  sqrt(M2/DataVector.n_elem);
	for (int_fast32_t n = 0; n < int_fast32_t(DataVector.n_elem); n++) {
		DataVector.at(n) = (DataVector.at(n) - mean)/popSD;
	}
	return Rcpp::wrap(DataVector);
}

static const R_CallMethodDef callMethods[] = {
	{"LocateLambdaCpp", (DL_FUNC) &LocateLambdaCpp, 3},
	{"LocateLambdaCpp_legacy", (DL_FUNC) &LocateLambdaCpp_legacy, 3},

	{"SkewVarCpp", (DL_FUNC) &SkewVarCpp, 2},
	{"SkewVar2Cpp", (DL_FUNC) &SkewVar2Cpp, 2},
	{"SkewAVarCpp", (DL_FUNC) &SkewAVarCpp, 2},

	{"colVarsCpp", (DL_FUNC) &colVarsCpp, 1},
	{"colSDsCpp", (DL_FUNC) &colSDsCpp, 1},
	{"colMeanSDCpp", (DL_FUNC) &colMeanSDCpp, 1},
	
	{"NZcolLog1pMeanSDCpp", (DL_FUNC) &NZcolLog1pMeanSDCpp, 2},
	{"NZcolLogMeanSDSkewCpp", (DL_FUNC) &NZcolLogMeanSDSkewCpp, 1},	
	{"NZcolMeansCpp", (DL_FUNC) &NZcolMeansCpp, 1},
	{"NZcolMeanSDCpp", (DL_FUNC) &NZcolMeanSDCpp, 1},
	{"WNZcolMeansCpp", (DL_FUNC) &WNZcolMeansCpp, 2},
	{"WcolMeansCpp", (DL_FUNC) &WcolMeansCpp, 2},
	{"popZScoreCPP", (DL_FUNC) &popZScoreCPP, 1},
	{"XPMCpp", (DL_FUNC) &XPMCpp, 1},
	{"tXPMCpp", (DL_FUNC) &tXPMCpp, 1},
	{"BatchEffectCpp", (DL_FUNC) &BatchEffectCpp, 4},
	{"getSlopeCpp", (DL_FUNC) &getSlopeCpp, 2},
	{"LinearRegressionCpp", (DL_FUNC) &LinearRegressionCpp, 2},
	{"LinearRegressionZeroCpp", (DL_FUNC) &LinearRegressionZeroCpp, 2},
	{"LinearRegressionFPCpp", (DL_FUNC) &LinearRegressionFPCpp, 4},
	{NULL, NULL, 0}
};

extern "C" {
	void R_init_Linnorm(DllInfo *info) { 
		R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	}
}
