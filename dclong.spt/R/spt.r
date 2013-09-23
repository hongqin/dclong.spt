#' @title Sequential Permutation Test
#'
#' @description Performs (multiple) sequential permutation tests.
#'
#' @param data a matrix, data frame or list.
#' If \code{data} is a matrix or data frame,
#' then each row of \code{data} corresponds to a gene/test.
#' If \code{data} is a list,
#' then each element of \code{data} corresponds to a gene/test.
#'
#' @param tsf0 a user-defined function for calculating the observed test statistics.
#'
#' @param tsf a user-defined function for permuting data and recalculating the test statistics. 
#' \code{tsf0} and \code{tsf} must be defined in the way such that bigger test statistics are more extreme.
#' These two functions must have the same signature. 
#' Both of them take the permuted data as the first argument. 
#' Extra arguments passing to these two functions are allowed.
#'
#' @param h the number of significant test statistics to hit before early termination.
#'
#' @param n the maximum number of permutations to use, including the observed one, i.e., 
#' at most n - 1 permutations will be sampled.  
#'
#' @param ... extra parameters (except the first argument \code{data}) 
#' to be passed to the user-defined test statistic functions \code{tsf0} and \code{tsf}.
#' 
#' @return an object of \code{\link[base]{class}} ``spt''. 
#' An object of class ``spt'' is a list containing the following components:
#' \item{p}{pvalue of sequential permutation test for each gene/test.}
#' \item{h}{see description in the Arguments section.}
#' \item{n}{see description in the Arguments section.}
#' Note that though \code{spt} performs sequential permutation test, 
#' you can get results for regular permutation test by setting \code{h >= n}.
#'
#' @seealso \code{\link{cxxwrapper}} which wraps C++ code for doing sequential permutation test.
#'
#' @export  
#' @examples
#' # download data
#' if(!file.exists('spt_data.rda')){
#'     download.file('http://dclong.github.io/media/spt/spt_data.rda', 'spt_data.rda')
#' }
#' load("spt_data.rda")
#' # center and standardize columns of marker
#' X = scale(t(marker))
#' # center and standardized rows of barley 
#' Y = t(scale(t(barley)))
#' #sequntial permutation test for no correlation between gene expression for each gene
#' #and the markers (it might take a while)
#' # for the purpose of illustration, only the first 100 genes are tested since it's time consuming.
#' spt.out.corr = spt(Y[1:100,], max_abs_corr0, max_abs_corr, 10, 1000, X=X)
#' 
#' # sequential permutation test for no difference between two treatment groups for each gene
#' # only the first 100 genes are tested since it's time consuming
#' spt.out.2gd = spt(leukemia[1:100,], abs_mean_diff0, abs_mean_diff, 10, 1000, n1=5)
#'

spt = function(data, tsf0, tsf, h, n, ...){
    if(is.matrix(data) || is.data.frame(data)){
        data_length = nrow(data)
        p = double(data_length)
        for(i in 1:data_length){
            p[i] = spt_one(data[i,], tsf0, tsf, h, n, ...) 
        }
        result = list(p=p, h=h, n=n)
        class(result) = "spt"
        return(result)
    }else if(is.list(data)){
        data_length = length(data)
        p = double(data_length)
        for(i in 1:data_length){
            p[i] = spt_one(data[[i]], tsf0, tsf, h, n, ...)
        }
        result = list(p = p, h = h, n = n)
        class(result) = "spt"
        return(result)
    }else{
        stop("data must be a matrix, data frame or list.")
    }
}

#‘ Calculate the sequencial p value for one gene/test.
spt_one = function(data, tsf0, tsf, h, n, ...){
    t0 = tsf0(data, ...)
    seq_count = 0
    for(i in 1:(n-1)){
        seq_count = seq_count + (tsf(data, ...) >= t0)
        if(seq_count >= h){
          return(seq_count / i)
        }
    }
    return((seq_count + 1)/n)
}
