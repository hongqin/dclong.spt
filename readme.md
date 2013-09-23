Permutation tests are commonly used. Usually it doesn't take much
time to do a permutation test, but it can be time-consuming if one have to
do many permutation tests (e.g., in gene data analysis) at the same time.
This package allows users to do sequential permutation tests which can save
lots of time. This package has implemented (both in R and C++ version) the
two sequential permutation tests discussed in the paper 
"Bancroft, Du and Nettleton (2013). 
Estimation of False Discovery Rate Using Sequential Permutation P­Values."  
It also allows users to do general sequential
permutation test by offering user-defined functions for calculating
observed and permuted test statistics. For people who are confident in
writing C++ code and willing to use C++ code to speed up computation, this
package tries all efforts to make this easy (thanks to Rcpp and C++
variadic template). Basically a user only has to offer C++ code for
calculating observed and permuted test statistics (similar to what he/she
would do using pure R code).
