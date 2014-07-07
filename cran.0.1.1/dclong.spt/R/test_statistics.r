#' @title User-defined Test Statistic Functions
#'
#' @description The functions end with "0" calculate the observed test statistics.
#' The corresponding versions that doesn't end with "0" permutes data and then recalculate the test statistics.
#' 
#' @export
#' @rdname test_statistics
#' @param data the data corresponding to a gene/test. 
#' In this case, 
#' it's a vector containing observations from two samples.
#' The first \code{n1} observations belong to the first sample,
#' and other observations correspond to the second sample.
#'
#' @param n1 size of the first sample.
#'
abs_mean_diff0 = function(data, n1){
    abs(mean(data[1:n1]) - mean(data[-(1:n1)]))
}
#' @export
#' @rdname test_statistics
abs_mean_diff = function(data, n1){
    abs_mean_diff0(sample(data, length(data), replace=F), n1)
}
#' @export
#' @rdname test_statistics
#' @param y the data corresponding to a gene/test, 
#' which is a vector in this case.
#' It must be centered and standardized.
#' 
#' @param X a matrix with columns being covariates. 
#' Each of the columns/covariates must be centered and standardized.
#' 
max_abs_corr0 = function(y, X){
    max(abs(t(y) %*% X))
}
#' @export
#' @rdname test_statistics
max_abs_corr = function(y, X){
    max_abs_corr0(sample(y, length(y), replace=F), X)
}

