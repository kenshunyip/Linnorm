/*
The MIT License (MIT) 
Copyright (c) 2016 by Shun Hang Yip <shunyip@bu.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#include <cstdlib>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <RcppArmadillo.h>
#include <vector>
#include <random>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

double gammaloglik(arma::vec vec3, double k, double adduplogvec, double addupvec) {	
	return (k-1)*adduplogvec - vec3.n_elem  * (k * (1 + log(addupvec/(vec3.n_elem*k)) ) + lgamma(k) );
}
// [[Rcpp::export]]
double gammaShape (arma::vec vec2) {
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


double SkewVarKurt	(const arma::mat& GeneExp, const double& lambda2) {
	arma::vec skewvec(GeneExp.n_rows);
	arma::vec SDevvec(GeneExp.n_rows);
	arma::vec kurtvec(GeneExp.n_rows);
	arma::vec meanvec(GeneExp.n_rows);

	double SumKurt = 0;
	double SumSkew = 0;
	double SumSD = 0;
	double SumMean = 0;
	double SumKurtMean = 0;
	double SumSkewMean = 0;
	double SumSDMean = 0;
	double SumMeanSq = 0;
	
	//Going to normalize mean to the scale from 1 to 2, need the range here to prevent looping all genes twice
	/*
	//double minMean = 0, maxMean = 0;
	//for (int i = 0; i < GeneExp.n_cols ; i++) {
	//	minMean += log(GeneExp.at(0,i) * lambda2 + 1);
	//	maxMean += log(GeneExp.at(GeneExp.n_rows-1,i) * lambda2 + 1);
	//}
	//double range = abs(maxMean - minMean);
	*/
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
		//Here, we can see that Variance is the second moment and Skewness is the third moment
		kurtvec.at(i) =  (GeneExp.n_cols * M4)/pow(M2,2) - 3;
		skewvec.at(i) =  (sqrt(GeneExp.n_cols) * M3) / pow(M2,1.5);
		SDevvec.at(i) =  sqrt(M2/(GeneExp.n_cols - 1));
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
	//y = Mx + C, here are the M and C
	double kurtM = (GeneExp.n_rows * SumKurtMean - SumKurt * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) );
	double kurtC = (SumKurt - kurtM * SumMean)/GeneExp.n_rows;
	double skewM = (GeneExp.n_rows * SumSkewMean - SumSkew * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) );
	double skewC = (SumSkew - skewM * SumMean)/GeneExp.n_rows;
	double SDevM = (GeneExp.n_rows * SumSDMean - SumSD * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) );
	
	//integral of the linear equation from 1 to length of geneexp
	double Skewzerointercept = -skewC/skewM;
	double Skewintegral = 0;
	if (Skewzerointercept > meanvec.at(0) && Skewzerointercept < meanvec.at(skewvec.n_elem-1)) {
		Skewintegral = abs((skewM * pow(Skewzerointercept,2))/2 + skewC * Skewzerointercept - (skewM * pow(meanvec.at(0),2))/2 - skewC * meanvec.at(0)) + abs((skewM * pow(meanvec.at(skewvec.n_elem-1),2))/2 + skewC * meanvec.at(skewvec.n_elem-1) - (skewM * pow(Skewzerointercept,2))/2 - skewC * Skewzerointercept);
	} else {
		Skewintegral = abs((skewM * pow(meanvec.at(skewvec.n_elem-1),2))/2 + skewC * meanvec.at(skewvec.n_elem-1) - (skewM * pow(meanvec.at(0),2))/2 - skewC * meanvec.at(0));
	}
	
	double Kurtzerointercept = -kurtC/kurtM;
	double Kurtintegral = 0;
	if (Kurtzerointercept > meanvec.at(0) && Kurtzerointercept < meanvec.at(kurtvec.n_elem-1)) {
		Kurtintegral = abs((kurtM * pow(Kurtzerointercept,2))/2 + kurtC * Kurtzerointercept - (kurtM * pow(meanvec.at(0),2))/2 - kurtC * meanvec.at(0)) + abs((kurtM * pow(meanvec.at(kurtvec.n_elem-1),2))/2 + kurtC * meanvec.at(kurtvec.n_elem-1) - (kurtM * pow(Kurtzerointercept,2))/2 - kurtC * Kurtzerointercept);
	} else {
		Kurtintegral = abs((kurtM * pow(meanvec.at(kurtvec.n_elem-1),2))/2 + kurtC * meanvec.at(kurtvec.n_elem-1) - (kurtM * pow(meanvec.at(0),2))/2 - kurtC * meanvec.at(0));
	}
	//cout << Skewintegral << " " <<  (Varintegral)/(SumMean/GeneExp.n_rows) << " " << Kurtintegral << endl;
	return (Skewintegral + Kurtintegral)/abs(meanvec.at(kurtvec.n_elem-1) - meanvec.at(0)) +  abs(SDevM);
}


double LocalSearch(const arma::mat& GeneExp, double minBound, double maxBound, double& smallest) {
	minBound = round(minBound);
	maxBound = round(maxBound);
	//cout << "here1 " << minBound<<   " " << maxBound << endl;
	double GminBound = minBound, GmaxBound = maxBound;
	double midBound = round((minBound + maxBound)/2);
	double om = abs(SkewVarKurt(GeneExp, midBound));
	double om2 = abs(SkewVarKurt(GeneExp,(midBound + 1)));
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
			minBound = round(smallestBound - (smallestBound - GminBound)/pow(2,index));
			maxBound = round(smallestBound + (GmaxBound - smallestBound)/pow(2,index));
			midBound = round((minBound + maxBound)/2);
			//cout << minBound << " " << maxBound << endl;
		}
		index++;
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

// [[Rcpp::export]]
double LocateLambda(const arma::mat& GeneExp, double search_exponent) {
	//Boundary of ILS are minBound and maxBound
	double leastMean = 0;
	int fivepercent = GeneExp.n_cols / 10;
	if (fivepercent < 1)
		fivepercent = 1;
	for (int j = 0; j < fivepercent; j++) {
		for (int i = 0; i < GeneExp.n_cols ; i ++) {
			leastMean += GeneExp.at(j,i);		
		}
	}
	leastMean = leastMean/(GeneExp.n_cols * fivepercent);
	double maxBound = 1/leastMean;
	double minBound = 1;
	//Find local min
	double localminIntegral;
	double localmin = LocalSearch(GeneExp, minBound, maxBound,localminIntegral);
	
	
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

double SkewVarKurt2	(const arma::mat& GeneExp, const double& lambda2) {
	arma::vec skewvec(GeneExp.n_rows);
	arma::vec SDevvec(GeneExp.n_rows);
	arma::vec kurtvec(GeneExp.n_rows);
	arma::vec meanvec(GeneExp.n_rows);

	double SumKurt = 0;
	double SumSkew = 0;
	double SumSD = 0;
	double SumMean = 0;
	double SumKurtMean = 0;
	double SumSkewMean = 0;
	double SumSDMean = 0;
	double SumMeanSq = 0;
	
	//Going to normalize mean to the scale from 1 to 2, need the range here to prevent looping all genes twice
	/*
	//double minMean = 0, maxMean = 0;
	//for (int i = 0; i < GeneExp.n_cols ; i++) {
	//	minMean += log(GeneExp.at(0,i) * lambda2 + 1);
	//	maxMean += log(GeneExp.at(GeneExp.n_rows-1,i) * lambda2 + 1);
	//}
	//double range = abs(maxMean - minMean);
	*/
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < GeneExp.n_rows; i ++) {
		double mean = 0;
		double M2 = 0;
		double M3 = 0;
		double M4 = 0;
		double delta, delta_n, delta_n2, term1;
		for (int n = 0; n < GeneExp.n_cols; n++) {
			delta = pow(GeneExp.at(i,n), lambda2) - mean;
			delta_n = delta / (n + 1);
			delta_n2 = delta_n * delta_n;
			term1 = delta * delta_n * n;
			mean = mean + delta_n;
			M4 = M4 + term1 * delta_n2 * (pow(n + 1,2) - 3 * (n + 1) + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
			M3 = M3 + term1 * delta_n * (n - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
		}
		//Here, we can see that Variance is the second moment and Skewness is the third moment
		kurtvec.at(i) =  (GeneExp.n_cols * M4)/pow(M2,2) - 3;
		skewvec.at(i) =  (sqrt(GeneExp.n_cols) * M3) / pow(M2,1.5);
		SDevvec.at(i) =  sqrt(M2/(GeneExp.n_cols - 1));
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
	//y = Mx + C, here are the M and C
	double kurtM = (GeneExp.n_rows * SumKurtMean - SumKurt * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) );
	double kurtC = (SumKurt - kurtM * SumMean)/GeneExp.n_rows;
	double skewM = (GeneExp.n_rows * SumSkewMean - SumSkew * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) );
	double skewC = (SumSkew - skewM * SumMean)/GeneExp.n_rows;
	double SDevM = (GeneExp.n_rows * SumSDMean - SumSD * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) );
	
	//integral of the linear equation from 1 to length of geneexp
	double Skewzerointercept = -skewC/skewM;
	double Skewintegral = 0;
	if (Skewzerointercept > meanvec.at(0) && Skewzerointercept < meanvec.at(skewvec.n_elem-1)) {
		Skewintegral = abs((skewM * pow(Skewzerointercept,2))/2 + skewC * Skewzerointercept - (skewM * pow(meanvec.at(0),2))/2 - skewC * meanvec.at(0)) + abs((skewM * pow(meanvec.at(skewvec.n_elem-1),2))/2 + skewC * meanvec.at(skewvec.n_elem-1) - (skewM * pow(Skewzerointercept,2))/2 - skewC * Skewzerointercept);
	} else {
		Skewintegral = abs((skewM * pow(meanvec.at(skewvec.n_elem-1),2))/2 + skewC * meanvec.at(skewvec.n_elem-1) - (skewM * pow(meanvec.at(0),2))/2 - skewC * meanvec.at(0));
	}
	
	double Kurtzerointercept = -kurtC/kurtM;
	double Kurtintegral = 0;
	if (Kurtzerointercept > meanvec.at(0) && Kurtzerointercept < meanvec.at(kurtvec.n_elem-1)) {
		Kurtintegral = abs((kurtM * pow(Kurtzerointercept,2))/2 + kurtC * Kurtzerointercept - (kurtM * pow(meanvec.at(0),2))/2 - kurtC * meanvec.at(0)) + abs((kurtM * pow(meanvec.at(kurtvec.n_elem-1),2))/2 + kurtC * meanvec.at(kurtvec.n_elem-1) - (kurtM * pow(Kurtzerointercept,2))/2 - kurtC * Kurtzerointercept);
	} else {
		Kurtintegral = abs((kurtM * pow(meanvec.at(kurtvec.n_elem-1),2))/2 + kurtC * meanvec.at(kurtvec.n_elem-1) - (kurtM * pow(meanvec.at(0),2))/2 - kurtC * meanvec.at(0));
	}
	//cout << Skewintegral << " " <<  (Varintegral)/(SumMean/GeneExp.n_rows) << " " << Kurtintegral << endl;
	return (Skewintegral + Kurtintegral)/abs(meanvec.at(kurtvec.n_elem-1) - meanvec.at(0)) +  abs(SDevM);
}


double LocalSearch2(const arma::mat& GeneExp, double minBound, double maxBound, double& smallest) {
	minBound = round(minBound * 10000)/10000;
	maxBound = round(maxBound* 10000)/10000;
	//cout << "here1 " << minBound<<   " " << maxBound << endl;
	double GminBound = minBound, GmaxBound = maxBound;
	double midBound = round((minBound + maxBound)/2);
	double om = abs(SkewVarKurt(GeneExp, midBound));
	double om2 = abs(SkewVarKurt(GeneExp,(midBound + 1)));
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
	midBound = round(((minBound + maxBound)/2)* 10000)/10000;
	int index = 0;
	while (smallestBound != midBound) {
		//cout << "here1 " << smallestBound<<   " " << midBound << endl;
		if (index > 0) {
			minBound = round((smallestBound - (smallestBound - GminBound)/pow(2,index))* 10000)/10000;
			maxBound = round((smallestBound + (GmaxBound - smallestBound)/pow(2,index))* 10000)/10000;
			midBound = round(((minBound + maxBound)/2)* 10000)/10000;
			//cout << minBound << " " << maxBound << endl;
		}
		index++;
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
			midBound = round(((minBound + maxBound)/2)* 10000)/10000;
		}
	}
	return midBound;
}

// [[Rcpp::export]]
double LocateLambda2( arma::mat GeneExp, double search_exponent) {
	//Boundary of ILS are minBound and maxBound
	double leastMean = 0;
	int fivepercent = GeneExp.n_cols / 10;
	if (fivepercent < 1)
		fivepercent = 1;
	for (int j = 0; j < fivepercent; j++) {
		for (int i = 0; i < GeneExp.n_cols ; i ++) {
			leastMean += GeneExp.at(j,i);		
		}
	}
	leastMean = leastMean/(GeneExp.n_cols * fivepercent);
	double maxBound = 1;
	double minBound = 0.00001;
	//Find local min
	double localminIntegral;
	double localmin = LocalSearch2(GeneExp, minBound, maxBound,localminIntegral);
	
	
	//LHS search
	int searchIndex = 1;
	//1.perturbing the current local minimum;
	double newminBound = localmin - search_exponent;
	double lastnewminBound = localmin;
	double newmin, newminIntegral;	
	while (newminBound >= minBound) {
		//cout << "LHS " << searchIndex <<  " " << newminBound <<endl;
		//2.applying local search after starting from the modified solution.
		newmin = LocalSearch2(GeneExp, newminBound, lastnewminBound, newminIntegral);
		//If new minimum is smaller than the local minimal, reset local minimal and start the process again.
		if (newminIntegral < localminIntegral) {
			localmin = newmin;
			localminIntegral = newminIntegral;
			lastnewminBound = newminBound;
			searchIndex++;
			//1.perturbing the current local minimum;
			if (newminBound != minBound) {
				newminBound = localmin -  search_exponent;
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
				newminBound = localmin -  search_exponent;
				if (newminBound < minBound) {
					newminBound = minBound;
				}
			} else {
				newminBound = minBound - 1;
			}
		}
	}
	//RHS search
	searchIndex = 1;
	//1.perturbing the current local maximum;
	double newmaxBound = localmin + search_exponent;
	double lastnewmaxBound = localmin;
	double newmax, newmaxIntegral;	
	while (newmaxBound <= maxBound) {
		//cout << "RHS " << searchIndex <<  " " << newmaxBound <<endl;
		//2.applying local search after starting from the modified solution.
		newmax = LocalSearch2(GeneExp, lastnewmaxBound, newmaxBound, newmaxIntegral);
		//If new maximum is smaller than the local maximal, reset local maximal and start the process again.
		if (newmaxIntegral < localminIntegral) {
			localmin = newmax;
			localminIntegral = newmaxIntegral;
			searchIndex++;
			lastnewmaxBound = newmaxBound;
			//1.perturbing the current local maximum;
			if (newmaxBound != maxBound) {
				newmaxBound = localmin +  search_exponent;
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
				newmaxBound = localmin + search_exponent;
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

