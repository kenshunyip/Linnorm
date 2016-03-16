/*
The MIT License (MIT) 
Copyright (c) <2016> <Shun Hang Yip>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <RcppArmadillo.h>
#include <vector>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

//Log likelihood function of gamma distribution
double gammaloglik(arma::vec vec3, double k, double adduplogvec, double addupvec) {	
	return (k-1)*adduplogvec - vec3.n_elem  * (k * (1 + log(addupvec/(vec3.n_elem*k)) ) + lgamma(k) );
}

//Maximum likelihood estimation for gamma distribution
// [[Rcpp::export]]
double gammaShape (arma::vec vec2) {
	//First, find the boundary to search from, based on slope.
	double minBound = 0, maxBound = 2;
	double midBound = (minBound + maxBound)/2 ;
	double adduplogvec = 0,addupvec = 0;
	for (int i = 0; i < vec2.n_elem; i++) {
		adduplogvec = adduplogvec + log(vec2.at(i));
		addupvec = addupvec + vec2.at(i);
	}
	double LL1 = gammaloglik(vec2,midBound,adduplogvec,addupvec);
	double LL2 = gammaloglik(vec2,maxBound,adduplogvec,addupvec);	
	while (LL1 < LL2) {
		minBound = midBound;
		maxBound = midBound * 2;
		midBound = (minBound + maxBound)/2 ;
		LL1 = gammaloglik(vec2,midBound,adduplogvec,addupvec);
		LL2 = gammaloglik(vec2,maxBound,adduplogvec,addupvec);
	}
	//Look for k using binary search.
	if (LL1 > LL2) {
		maxBound = midBound;
		midBound = (minBound + maxBound)/2 ;
		while (midBound != minBound) {
			LL1 = gammaloglik(vec2,midBound,adduplogvec,addupvec);
			LL2 = gammaloglik(vec2,maxBound,adduplogvec,addupvec);
			if (LL1 > LL2) {
				maxBound = midBound;
			} else {
				minBound = midBound;
			}
			midBound = (minBound + maxBound)/2 ;
		}
	} else {
		//Safety here. If likelihood can go up to infinity, we find the smallest k where the program can't tell the difference between LL1 and LL2 anymore.
		while(LL1 == LL2) {
			minBound = midBound/2;
			maxBound = midBound;
			midBound = (minBound + maxBound)/2 ;
			LL1 = gammaloglik(vec2,midBound,adduplogvec,addupvec);
			LL2 = gammaloglik(vec2,maxBound,adduplogvec,addupvec);
		}
		while(midBound != minBound && midBound != maxBound) {
			if (LL1 != LL2) {
				minBound = midBound;
			} else {
				maxBound = midBound;
			}
			midBound = (minBound + maxBound)/2 ;
			LL1 = gammaloglik(vec2,midBound,adduplogvec,addupvec);
			LL2 = gammaloglik(vec2,maxBound,adduplogvec,addupvec);
		}		
	}
	
	return midBound;
}

//This function calculates F(lambda) (see article) based on the expression matrix and lambda.
double SkewVarKurt	(const arma::mat& GeneExp, const double& lambda2) {
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::vec skewvec(GeneExp.n_rows);
	arma::vec SDevvec(GeneExp.n_rows);
	arma::vec kurtvec(GeneExp.n_rows);
	arma::vec meanvec(GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	double SumKurt = 0;
	double SumSkew = 0;
	double SumSD = 0;
	double SumMean = 0;
	double SumKurtMean = 0;
	double SumSkewMean = 0;
	double SumSDMean = 0;
	double SumMeanSq = 0;
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < GeneExp.n_rows; i ++) {
		double mean = 0;
		double M2 = 0;
		double M3 = 0;
		double M4 = 0;
		double delta, delta_n, delta_n2, term1;
		for (int n = 0; n < GeneExp.n_cols; n++) {
			delta = log1p(GeneExp.at(i,n) * lambda2) - mean;
			delta_n = delta / (n + 1);
			delta_n2 = delta_n * delta_n;
			term1 = delta * delta_n * n;
			mean = mean + delta_n;
			M4 = M4 + term1 * delta_n2 * (pow(n + 1,2) - 3 * (n + 1) + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
			M3 = M3 + term1 * delta_n * (n - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
		}
		//Here, calculate kurtosis, skewness and SD using 4th, 3rd and 2nd moments.
		kurtvec.at(i) =  (GeneExp.n_cols * M4)/pow(M2,2) - 3;
		skewvec.at(i) =  (sqrt(GeneExp.n_cols) * M3) / pow(M2,1.5);
		SDevvec.at(i) =  sqrt(M2/(GeneExp.n_cols - 1));
		
		//Linear regression
		//meanvec.at(i) = 1 + (mean - minMean)/range;
		//meanvec.at(i) = 1 + i;
		meanvec.at(i) = mean;
		SumKurt += kurtvec.at(i);
		SumSkew += skewvec.at(i);
		SumSD += SDevvec.at(i);
		SumMean += meanvec.at(i);
		SumKurtMean += kurtvec.at(i) * meanvec.at(i);
		SumSkewMean += skewvec.at(i) * meanvec.at(i);
		SumSDMean += SDevvec.at(i) * meanvec.at(i);
		SumMeanSq += pow(meanvec.at(i),2);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	double kurtM = (GeneExp.n_rows * SumKurtMean - SumKurt * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) ) ;
	double kurtC = (SumKurt - kurtM * SumMean)/GeneExp.n_rows;
	double skewM = (GeneExp.n_rows * SumSkewMean - SumSkew * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) ) ;
	double skewC = (SumSkew - skewM * SumMean)/GeneExp.n_rows;
	double SDevM = (GeneExp.n_rows * SumSDMean - SumSD * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) );
	
	//integral of the linear equation of skewness
	double Skewzerointercept = -skewC/skewM;
	double Skewintegral = 0;
	skewM = skewM/2;
	if (Skewzerointercept > meanvec.at(0) && Skewzerointercept < meanvec.at(skewvec.n_elem-1)) {
		Skewintegral = (abs(skewM * (meanvec.at(skewvec.n_elem-1) + Skewzerointercept) + skewC) * (meanvec.at(skewvec.n_elem-1) - Skewzerointercept) + abs(skewM * (meanvec.at(0) + Skewzerointercept) + skewC)* (Skewzerointercept - meanvec.at(0)))/(meanvec.at(kurtvec.n_elem-1) - meanvec.at(0));
	} else {
		Skewintegral = abs(skewM * (meanvec.at(skewvec.n_elem-1) + meanvec.at(0)) + skewC);
	}
	//integral of the linear equation of kurtosis
	double Kurtzerointercept = -kurtC/kurtM;
	double Kurtintegral = 0;
	kurtM = kurtM/2;
	if (Kurtzerointercept > meanvec.at(0) && Kurtzerointercept < meanvec.at(kurtvec.n_elem-1)) {
		Kurtintegral = (abs(kurtM * (meanvec.at(kurtvec.n_elem-1) + Kurtzerointercept) + kurtC) * (meanvec.at(kurtvec.n_elem-1) - Kurtzerointercept) + abs(kurtM * (meanvec.at(0) + Kurtzerointercept) + kurtC)* (Kurtzerointercept - meanvec.at(0)) )/(meanvec.at(kurtvec.n_elem-1) - meanvec.at(0));
	} else {
		Kurtintegral = abs(kurtM * (meanvec.at(kurtvec.n_elem-1) + meanvec.at(0)) + kurtC);
	}
	return Skewintegral + Kurtintegral +  abs(SDevM);
}

//Given a range of lambda, this function finds lambda that minimizes F(lambda) (see article) based on the expression matrix by using binary search.
double LocalSearch(const arma::mat& GeneExp, double minBound, double maxBound, double& smallest) {
	minBound = round(minBound);
	maxBound = round(maxBound);
	//cout << "here1 " << minBound<<   " " << maxBound << endl;
	double GminBound = minBound, GmaxBound = maxBound;
	double midBound = round((minBound + maxBound)/2);
	double om = abs(SkewVarKurt(GeneExp, midBound));
	double om2 = abs(SkewVarKurt(GeneExp,(midBound + 1)));
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
	int index = 0;
	while (smallestBound != midBound) {
		//cout << "here1 " << smallestBound<<   " " << midBound << endl;
		if (index > 0) {
			//Reset boundary to center at smallestBound
			minBound = round(smallestBound - (smallestBound - GminBound)/pow(2,index));
			maxBound = round(smallestBound + (GmaxBound - smallestBound)/pow(2,index));
			midBound = round((minBound + maxBound)/2);
			//cout << minBound << " " << maxBound << endl;
		}
		index++;
		//Binary search for smallest F(lambda)
		while (minBound != midBound && maxBound != midBound) {
			om = abs(SkewVarKurt(GeneExp, midBound));
			om2 = abs(SkewVarKurt(GeneExp,midBound + 1));
			if (om > om2) {
				minBound = midBound;
				if (om2 < smallest) {
					smallest = om2;
					smallestBound = midBound + 1;
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
	return midBound;
}

//Using the LocalSearch function, performs iterated local search to find minimal lambda.
// [[Rcpp::export]]
double LocateLambda(const arma::mat& GeneExp, double search_exponent) {
	//Boundary of ILS are minBound and maxBound
	//Here, we define the range of lambda based on the dataset.
	double leastMean = 0;
	int fivepercent = GeneExp.n_cols / 20;
	if (fivepercent < 1)
		fivepercent = 1;
	int numnonzero = 0;
	for (int j = 0; j < fivepercent; j++) {
		for (int i = 0; i < GeneExp.n_cols ; i ++) {
			if ( GeneExp.at(j,i) > 0) {
				leastMean += GeneExp.at(j,i);
				numnonzero++;
			}				
		}
	}
	leastMean = leastMean/numnonzero;
	double maxBound = 1/leastMean;
	double minBound = 1;
	
	//Find local minima
	double localminIntegral;
	double localmin = LocalSearch(GeneExp, minBound, maxBound,localminIntegral);
	
	//First, iterated local search starting on the left hand side of the local minima.  
	//LHS search
	int searchIndex = 1;
	//1.perturbing the current local minimum;
	double newminBound = localmin - 10 * search_exponent;
	double lastnewminBound = localmin;
	double newmin, newminIntegral;	
	while (newminBound >= minBound) {
		//cout << "LHS " << searchIndex <<  " " << newminBound <<endl;
		//2.applying local search after starting from the modified solution.
		newmin = LocalSearch(GeneExp, newminBound, lastnewminBound, newminIntegral);
		//If new minimum is smaller than the local minimal, reset local minimal and start the process again.
		if (newminIntegral < localminIntegral) {
			localmin = newmin;
			localminIntegral = newminIntegral;
			lastnewminBound = newminBound;
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
			lastnewminBound = newminBound;
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
	searchIndex = 1;
	//1.perturbing the current local maximum;
	double newmaxBound = localmin + 10 * search_exponent;
	double lastnewmaxBound = localmin;
	double newmax, newmaxIntegral;	
	while (newmaxBound <= maxBound) {
		//cout << "RHS " << searchIndex <<  " " << newmaxBound <<endl;
		//2.applying local search after starting from the modified solution.
		newmax = LocalSearch(GeneExp, lastnewmaxBound, newmaxBound, newmaxIntegral);
		//If new maximum is smaller than the local maximal, reset local maximal and start the process again.
		if (newmaxIntegral < localminIntegral) {
			localmin = newmax;
			localminIntegral = newmaxIntegral;
			searchIndex++;
			lastnewmaxBound = newmaxBound;
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
			lastnewmaxBound = newmaxBound;
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
	return localmin;
}

