/*
 * spt.h
 *
 *  Created on: Feb 6, 2012
 *      Author: adu
 */

#ifndef SPT_H_
#define SPT_H_


//void rpt_corr(double*x, int*covNum, int*covLength, double*y, int*resNum, int*n,
//		double*pvalues, double*obsMax, int*maxIndex);
//
//void rpt1_corr(double*x, int*covNum, int*covLength, double*xsd, double*y, double*ysd, int*n,
//		double*pvalue, double*obsMax, int*maxIndex);

void spt_corr(double*x, int*xrow, int*xcol, double*y, int*yrow, int*h, int*n,
		double*pvalue, double*obsMax, int*maxIndex);

void spt1_corr(double*x, int*xrow, int*xcol, double*xsd, double*y, double*ysd, int*h,
		int*n, double*pvalue, double*obsMax, int*maxIndex);

void maxAbsCor(double*x, int*covNum, int*xlength, double*xsd, double*y,
		double*ysd, double*max, int*maxIndex);

double cor(double*x, double*xsd, double*y, double*ysd, int*n);

double scale(double*v, int*n);

void shuffle(double*x, int*n);
/**
 * Randomly shuffle vector x to generate a random permutation.
 * Note that this function mutate the original array.
 * @param x a pointer to a double array to be shuffled.
 * @param n a pointer to the length of vector x.
 */

void shuffle(double*x, int*n) {
	int index, i;
	double temp;
	for (i = *n - 1; i > 0; i--) {
		//generate an integer in 0 and n-1 (inclusive)
		index = (int) runif(0, i + 1);
		//no need, but just in case
		if (index > i) {
			index = i;
		}
		//swap x[i] and x[index]
		temp = x[i];
		x[i] = x[index];
		x[index] = temp;
	}
}

//void rpt_corr(double*x, int*covNum, int*covLength, double*y, int*resNum, int*n,
//		double*pvalues, double*obsMax, int*maxIndex) {
//	//scale all vectors
//	double* xsd = (double*) malloc(*covNum * sizeof(double));
//	double* ysd = (double*) malloc(*resNum * sizeof(double));
//	int i;
//	for (i = 0; i < *covNum; i++) {
//		//center the independent variables and calculate their proportional standard deviation for future use
//		xsd[i] = scale(x + i * (*covLength), covLength);
//	}
//	//get the RNG state in R
//	GetRNGstate();
//	for (i = 0; i < *resNum; i++) {
//		//center the response vector and calculate its proportional standard deviation for future use
//		ysd[i] = scale(y + i * (*covLength), covLength);
//		//perform regular permutation test for this response vector
//		rpt1_corr(x, covNum, covLength, xsd, y + i * (*covLength), ysd + i, n,
//						pvalues + i, obsMax + i, maxIndex + i);
//	}
//	//write the RNG state back to R
//	PutRNGstate();
//	free(xsd);
//	free(ysd);
//}

/**
 * This function performs a sequential permutation test.
 * It is meant to be called from R, so its parameters are all pointers.
 * @param x a pointer to a long array containing all covariates.
 * @param covNum a pointer to the number of rows of x.
 * @param covLength a pointer to the number of columns of x.
 * @param y a pointer to another long array containing all response variables.
 * @param resNum a pointer to the number of rows of y. Notice that it's unnecessary to specify the number of columns for
 * y since it's the same as *covNum.
 * @param h a pointer to the stop criteria.
 * @param n a pointer to the simulation size.
 * @param pvalue a pointer to the returned pvalues whose length is *yrow (i.e. the same as the number of ROWs of y).
 * @param obsMax a pointer to the returned maximum absolute correlations.
 * @param maxIndex a pointer to the returned index of the maximum absolute correlations.
 */
void spt_corr(double *x, int* covNum, int* covLength, double* y, int* resNum, int* h,
		int* n, double* pvalues, double*obsMax, int*maxIndex) {
	//scale all vectors
	double* xsd = (double*) malloc(*covNum * sizeof(double));
	double* ysd = (double*) malloc(*resNum * sizeof(double));
	int i;
	for (i = 0; i < *covNum; i++) {
		//center the independent variables and calculate their proportional standard deviation for future use
		xsd[i] = scale(x + i * (*covLength), covLength);
	}
	//get the RNG state in R
	GetRNGstate();
	for (i = 0; i < *resNum; i++) {
		//center the response vector and calculate its proportional standard deviation for future use
		ysd[i] = scale(y + i * (*covLength), covLength);
		//perform sequential permutation test for this response vector
		spt1_corr(x, covNum, covLength, xsd, y + i * (*covLength), ysd + i, h, n,
						pvalues + i, obsMax + i, maxIndex + i);
	}
	//write the RNG state back to R
	PutRNGstate();
	free(xsd);
	free(ysd);
}

//void rpt1_corr(double*x, int*covNum, int*covLength, double*xsd, double*y, double*ysd, int*n,
//		double*pvalue, double*obsMax, int*maxIndex) {
//	double permCor;
//	int i, j;
//	//initial *pvalue
//	*pvalue = 0;
//	//get the observed maximum absolute correlation and its index
//	maxAbsCor(x, covNum, covLength, xsd, y, ysd, obsMax, maxIndex);
//	for (i = 1; i <= *n; i++) {
//		//randomly shuffle y
//		rshuffle(y, covLength,covLength);
//		for (j = 0; j < *covNum; j++) {
//			permCor = cor(x + j * (*covLength), xsd + j, y, ysd, covLength);
//			if (permCor < 0) {
//				permCor *= -1;
//			}
//			if (permCor >= *obsMax) {
//				(*pvalue)++;
//				break;
//			}
//		}
//	}
//	//increase *pvalue by 1 to include to observed one
//	(*pvalue)++;
//	*pvalue /= (*n + 1);
//}
/**
 * Performs a sequential permutation test for a correlation between a single y (response variable) and a bunch of x's (covariates).
 * Note that these covariates are passed together in a long array.
 * @param x a pointer to a the array that contains all covariates.
 * @param covNum a pointer to the number of rows of x.
 * @param covLength a pointer to the number of columns of x.
 * @param xsd a pointer to the proportional standard deviations of x.
 * @param y a pointer to a double array.
 * @param ysd a pointer to the proportional standard deviations of y.
 * @param h a pointer to the stop criteria.
 * @param n a pointer to the simulation size.
 * @param pvalue a pointer to store the sequential permutation test pvalue.
 * @param maxIndex a pointer to store the index of covariate at which the maximum correlation is achieved.
 */
void spt1_corr(double *x, int* covNum, int* covLength, double*xsd, double *y,
		double *ysd, int*h, int*n, double*pvalue, double*obsMax, int*maxIndex) {
	double permCor;
	int i, j;
	//initial *pvalue
	*pvalue = 0;
	//get the observed maximum absolute correlation and its index
	maxAbsCor(x, covNum, covLength, xsd, y, ysd, obsMax, maxIndex);
	for (i = 1; i <= *n; i++) {
		//randomly shuffle y
		shuffle(y, covLength);
		for (j = 0; j < *covNum; j++) {
			permCor = cor(x + j * (*covLength), xsd + j, y, ysd, covLength);
			if (permCor < 0) {
				permCor *= -1;
			}
			if (permCor >= *obsMax) {
				(*pvalue)++;
				//check whether h correlations has exceeded obsMax
				if (*pvalue == *h) {
					//return Monte Carlo pvalue
					*pvalue /= i;
					return;
				}
				break;
			}
		}
	}
	//increase *pvalue by 1 to include to observed one
	(*pvalue)++;
	*pvalue /= (*n + 1);
}

/**
 * Calculates the maximum absolute correlation between a vector y and a bunch of covariates x's.
 * Notice that these covariates are passed together as a long array.
 * @param x a pointer to a the array containing all covariates.
 * @param covNum a pointer to the number of covariates in x.
 * @param covLength a pointer to the length of each covariate.
 * @param xsd a pointer to the array containing standard deviations of x's.
 * @param y a pointer to a 1-D double array.
 * @param ysd a pointer to the standard deviation of y.
 * @param max a pointer to store the maximum absolute correlation.
 * @param maxIndex a pointer to store the index of the covariate at which the maximum correlation is achieved.
 */
void maxAbsCor(double*x, int*covNum, int*covLength, double*xsd, double*y,
		double*ysd, double*max, int*maxIndex) {
	double current;
	int i;
	//initialize *max and *maxIndex
	*max = 0;
	*maxIndex = -1;
	for (i = 0; i < *covNum; i++) {
		current = cor(x + i * (*covLength), xsd + i, y, ysd, covLength);
		if (current < 0) {
			current *= -1;
		}
		if (current > *max) {
			*max = current;
			*maxIndex = i;
		}
	}
}

/**
 * Calculates the correlation between vector x and y.
 * This function only works on vectors that are centered, i.e., it should
 * only be used after function scale has been applied to these vectors.
 * @param x a pointer to a vector.
 * @param xsd a pointer to the proportional standard deviation of vector x.
 * @param y a pointer to another vector.
 * @param ysd a pointer to the proportional standard deviation of vector y.
 * @param n a pointer to the length of the two vectors.
 * @return the correlation between the two vectors.
 */
double cor(double*x, double*xsd, double*y, double*ysd, int*n) {
	double correlation = 0;
	int i;
	for (i = 0; i < *n; i++) {
		correlation += x[i] * y[i];
	}
	return (correlation / *xsd / *ysd);
}

/**
 * Center and a vector calculate its (proportional) standard deviation for later use.
 * @param v a pointer to a vector.
 * @param n a pointer to the length of the vector.
 * @return the proportional standard deviation of the vector.
 */
double scale(double*v, int*n) {
	//center the vector
	double ave = 0;
	double sd = 0;
	int i;
	for (i = 0; i < *n; i++) {
		ave += v[i];
	}
	ave /= *n;
	for (i = 0; i < *n; i++) {
		v[i] -= ave;
		sd += v[i] * v[i];
	}
	return sqrt(sd);
}


#endif /* SPT_H_ */
