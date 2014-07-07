#' @title Sequential Permutation Test
#'
#' @description Performs (multiple) sequential permutation tests for no correlations between response variables and covariates.
#'
#' @param x a matrix/vector with columns containing covariates.
#' @param y a matrix/vector with columns containing response variables.
#' @param h the number of significant test statistics to hit before early termination.
#' @param n the maximum number of permutations to use, including the observed one, i.e., 
#' at most n - 1 permutations will be sampled.  
#'
#' @details For each response varaible \code{y0} (a column in \code{y}), 
#' a sequential permutation test is done for H0: there's correlation between \code{y0} 
#' and covariates (columns in \code{x}) VS Ha: no correlation. 
#' The maximum absolute correlation is the test statistic used.  
#'
#' @return an object of \code{\link[base]{class}} ``spt''. 
#' An object of class ``spt'' is a list containing at least the following components:
#' \item{p}{pvalue of sequential permutation test for each gene/test.}
#' \item{h}{see description in the Arguments section.}
#' \item{n}{see description in the Arguments section.}
#' The object return by function \code{spt.corr} also contains the following component(s):
#' \item{max.ac}{the maximum absolute correlation (which is the test statistic) for each gene.}
#' \item{max.index}{the index of the covariate where the maximum absolute correlations is achieved for each gene/test.}
#' Note that though \code{spt.corr} performs sequential permutation test, you can
#' get results for regular permutation test by setting \code{h >= n}.
#'
#' @seealso \code{\link{spt.mean}} sequential permutation tests for 
#' no difference between two treatment groups.
#' @export  
#' @examples
#' \dontrun{
#' #load data
#' data(marker)
#' data(barley)
#' #sequntial permutation test for no correlation between gene expression
#' #and the markers (it might take a while)
#' spt.corr(t(marker),t(barley),10,1000)-> spt.corr.out
#' head(spt.corr.out)
#'}
spt.corr = function(x,y,h,n){
  covLength = nrow(x)
  if(covLength!=nrow(y)){
      stop('matrix x and y must have the same number of rows.')
  }
  covNum = ncol(x)
  resNum = ncol(y)
  result = .C("spt_corr",as.double(x),as.integer(covNum),as.integer(covLength),as.double(y),
     as.integer(resNum),as.integer(h),as.integer(n-1),
     double(resNum),double(resNum),integer(resNum),PACKAGE="dclong.spt")
  list(p=result[[8]],max.ac=result[[9]],max.index=result[[10]]+1,h=h,n=n)->result
  class(result) = "spt"
  return(result)
}

# rpt.corr = function(x,y,n){
#     covLength = nrow(x)
#     if(covLength!=nrow(y)){
#         stop('matrix x and y must have the same number of rows.')
#     }
#     covNum = ncol(x)
#     resNum = ncol(y)
#     result = .C("rpt_corr",as.double(x),as.integer(covNum),as.integer(covLength),as.double(y),
#         as.integer(resNum),as.integer(n),
#         double(resNum),double(resNum),integer(resNum))
#     list(pvalues=result[[7]],maxAbsCor=result[[8]],maxIndex=result[[9]]+1,n=n)
# }






