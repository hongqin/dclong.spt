#' @title Sequential Permutation Test 
#'
#' @description Performs (multiple) sequential permutation tests for 
#' no difference between two treatment groups.
#'
#' @param data a matrix with each row being a data set for a test (gene). 
#' @param size1 the size of the first group.
#' @inheritParams spt.corr 
#'
#' @details For each data set (a row in \code{data}) containing two groups, 
#' a sequential permutation test is done for H0: there's difference between 
#' the two groups VS no difference. The absolute difference between the means of 
#' the two groups is the test statistic used. 
#'
#' @return an object of \code{\link[base]{class}} ``spt''. 
#' An object of class ``spt'' is a list containing at least the following components:
#' \item{p}{pvalues of sequential permutation tests.}
#' \item{h}{see description in the Arguments section.}
#' \item{n}{see description in the Arguments section.}
#' The object return by function \code{spt.corr} also contains the following component(s):
#' \item{n.perms}{number of sequential permutations sampled for each test.}
#' Note that though \code{spt.mean} performs sequential permutation tests, 
#' you can get results for regular permutation tests by setting \code{h >= n}.
#'
#' @seealso \code{\link{spt.corr}} for sequential permutations test for no correlation
#' between a response variable and a bunch of covariates.
#' @export
#' @examples
#' \dontrun{
#'    #load data
#'    data(leukemia)
#'    spt.mean(leukemia,5,10,1000) -> spt.mean.out
#'    head(spt.mean.out)
#' }
spt.mean = function(data,size1,h,n){
  if(is.matrix(data)){
    numberOfGenes = nrow(data)
    numberOfObs = ncol(data)
    result = .C("spt_mean",as.double(t(data)),as.integer(numberOfGenes),as.integer(numberOfObs),
       as.integer(size1),as.integer(h),as.integer(n-1),
       double(numberOfGenes),double(numberOfGenes),integer(numberOfGenes),PACKAGE="dclong.spt")
    list(p=result[[8]],n.perms = result[[9]],h=h,n=n) -> result
    class(result) = "spt"
    return(result)
  }
  stop("argument data must be a matrix.")
}

