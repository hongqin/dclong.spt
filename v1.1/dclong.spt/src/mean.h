/*
 * mean.h
 *
 *  Created on: Feb 6, 2012
 *      Author: adu
 */

#ifndef MEAN_H_
#define MEAN_H_


void spt_mean(double*dataForGenes, int*numberOfGenes, int*numberOfObs,
		int*numberOfObsInFirstGroup, int*seqPermStopCriteria,
		int*numberOfPerms, double*rpvalues, double*spvalues,
		int*numberOfSeqPerms);

void spt1_mean(double*dataForGene, int n,
		int size1,int size2,int n1,int n2, int seqPermStopCriteria,
		int numberOfPerms, double*rpvalue, double*spvalue, int*numberOfSeqPerms);

double test_statistic(double*x,double sum, int n1, int n2);
/**
 * This function performs a sequential permutation test for no
 * difference between two groups for each gene using the difference between the means of
 * the two groups as the test statistic.
 * It is meant to be called from R, so its parameters are all pointers.
 * @param dataForGenes a pointer to a long array containing (two-group data) for all genes.
 * @param numberOfGenes a pointer to the number of genes.
 * @param numberOfObs a pointer to the number of observations in each gene.
 * @param numberOfObsInFirstGroup a pointer to the number of observations in the first group in each gene.
 * @param seqPermStopCriteria the number of significant test statistics to hit before early termination.
 * @param numberOfPerms a pointer to the number of (regular) permutations to be done.
 * @param rpvalues a pointer to the regular permutation pvalues.
 * @param spvalues a pointer to the sequential permutation test pvalues.
 * @param numberOfSeqPerm a pointer to the total number of sequential permutations.
 */
void spt_mean(double*dataForGenes, int*numberOfGenes, int*numberOfObs,
		int*numberOfObsInFirstGroup, int*seqPermStopCriteria,
		int*numberOfPerms, double*rpvalues, double*spvalues,
		int*numberOfSeqPerms) {
	int i;
	//-----------------------------------------------------
	int n = *numberOfObs;
	int size1 = *numberOfObsInFirstGroup;
	int size2 = n - size1;
	int n1, n2;
	if(size1<size2){
		n1 = size1;
		n2 = size2;
	}else{
		n1 = size2;
		n2 = size1;
	}
	//-----------------------------------------------------
	//get the RNG state in R
	GetRNGstate();
	//loop through genes
	for (i = 0; i < *numberOfGenes; ++i) {
		spt1_mean(dataForGenes + i * n, n,size1,size2,n1,n2,
				*seqPermStopCriteria, *numberOfPerms,
				rpvalues + i, spvalues + i, numberOfSeqPerms + i);
	}
	//write the RNG state back to R
	PutRNGstate();
}
/**
 * This function performs a sequential permutation test for no
 * difference between two groups for a single gene using the absolute difference between the means
 * of the two groups as the test statistic.
 * @param dataForGene a pointer to an array containing (two-group data) for the gene.
 * @param n the number of observations in the gene.
 * @param size1 the number of observations in the first group.
 * @param size2 the number of observations in the second group.
 * @param n1 the smaller one between size1 and size2.
 * @param n2 the bigger one between size1 and size2.
 * @param seqPermStopCriteria the criteria for terminate sequential permutation.
 * @param numberOfPerms the maximum number of permutations.
 * @param rpvalue a pointer to the regular permutation pvalue.
 * @param spvalue a pointer to the sequential permutation test pvalue.
 * @param numberOfSeqPerm a pointer to the total number of sequential permutations.
 */
void spt1_mean(double*dataForGene, int n,
		int size1,int size2,int n1,int n2, int seqPermStopCriteria,
		int numberOfPerms, double*rpvalue, double*spvalue, int*numberOfSeqPerms) {
	//initial variables;
	int i, j;
	double sum1 = 0;
	double sum2 = 0;
	double sum, absMeanDiff, tempAbsMeanDiff;
	int regularCount = 0;
	int sequentialCount = 0;
	*numberOfSeqPerms = 0;
	//sum of the first group
	for (i = 0; i < size1; ++i) {
		sum1 += dataForGene[i];
	}
	//sum of the second group
	for (i = size1; i < n; ++i) {
		sum2 += dataForGene[i];
	}
	//sum of the two groups
	sum = sum1 + sum2;
	//absolute difference between the means of the two groups
	absMeanDiff = sum1 / size1 - sum2 / size2;
	if(absMeanDiff < 0){
		absMeanDiff = - absMeanDiff;
	}
	//regular and sequential permutation tests
	for (i = 0; i < numberOfPerms; ++i) {
		tempAbsMeanDiff = test_statistic(dataForGene,sum,n1,n2);
		if (tempAbsMeanDiff >= absMeanDiff) {
			//increase count for the regular permutation test
			regularCount++;
			//increase count for the sequential permutation test
			if (sequentialCount < seqPermStopCriteria) {
				sequentialCount++;
				(*numberOfSeqPerms)++;
			}
		} else {
			if (sequentialCount < seqPermStopCriteria) {
				(*numberOfSeqPerms)++;
			}
		}
	}
	//regular permutation test pvalue
	regularCount++;
	numberOfPerms++;
	(*rpvalue) = regularCount / (double)numberOfPerms;
	//sequential permutation test pvalue
	if (sequentialCount == seqPermStopCriteria) {
		(*spvalue) = sequentialCount / (double)*numberOfSeqPerms;
	} else {
		sequentialCount++;
		(*spvalue) = sequentialCount / (double)numberOfPerms;
	}
}

double test_statistic(double*x,double sum, int n1, int n2){
	double sum1 = 0;
	int index, i;
	double temp;
	for(i=n1+n2-1; i>=n2; --i){
		index = (int)runif(0,i+1);
		if(index>i){
			index = i;
		}
		//swap x[i] and x[index]
		temp = x[i];
		x[i] = x[index];
		x[index] = temp;
		//sum the selected element
		sum1 += x[i];
	}
	temp = sum1/n1 - (sum-sum1)/n2;
	if(temp<0){
		temp = -temp;
	}
	return temp;
}


#endif /* MEAN_H_ */
