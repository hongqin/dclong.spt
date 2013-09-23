
#' @title C++ Code for Sequential Permutation Test
#'
#' @description The function \code{cxxwrapper} allows user 
#' to run C++ code for sequential permutation test easily.
#' Similar to what you can do with the R code in this package,
#' you can write C++ functions to permute data and calculate test statistics,
#' and then use the function \code{cxxwrapper} to compile your C++ code and 
#' generate an R function calling the compiled code. 
#' To solve the problem of passing data between R and C++, 
#' the underlying C++ code limits the data structure that you can use. 
#' \code{arma::mat} (a matrix with double elements) is currently the only supported data structure 
#' for storing data corresponding to all genes. 
#' Each row of the matrix contains data corresponding to one gene.  
#'
#' @param tsf_code the C++ code for user-defined test statistic functions. 
#' Similar to the R code, 
#' you have to write 2 C++ functions. 
#' One of them calculate the observed test statistic (based on original data),
#' the other permutes data and calculates the test statistic again. 
#' 
#' @param tsf0 the name of the C++ function which calculates the observed test statistic 
#' (based on original data). 
#' This function must take a (first) argument of type \code{arma::rowvec}. 
#' Extra parameters passing to this function are allowed. 
#' See \code{epList} for more on passing extra parameters to \code{tsf0}.
#' For convenience of user,
#' "abs_mean_diff0" and "max_abs_corr0" are defined in "rcpp/test_statistics.cpp"
#' under the home directory of this package. 
#'
#' @param tsf the name of the C++ function which permutes data and then calculate the test statistic again.  
#' \code{tsf0} and \code{tsf} must be defined in the way such that bigger test statistics are more extreme. 
#' This function must take a (first) argument of type \code{arma::rowvec}.
#' Extra parameters passing to this function are allowed.
#' See \code{epList} for more on passing extra parameters to \code{tsf}.
#' For convenience of user, 
#' "abs_mean_diff" and "max_abs_corr" are defined in "rcpp/test_statistics.cpp"
#' under the home directory of this package.
#' Using \code{abs_mean_diff0} and \code{abs_mean_diff}, 
#' one can do sequential permutation test for no loation difference between two samples. 
#' Using \code{max_abs_corr0} and \code{max_abs_corr},
#' one can do sequential permutation test for no correlation between a reponse variable and 
#' a bunch of covariates. 
#' 
#' @param epList a string matrix (or string vector if only 1 extra parameter) 
#' containing information of extra parameters 
#' to be passed to the user-defined functions \code{tsf0} and \code{tsf}.
#' Each row of \code{epList} containing information of an extra parameter. 
#' The first column of \code{epList} contains names of extra parameters. 
#' The names of extra parameters must be valid C++ and R variable names. 
#' The names cannot be \code{data}, \code{h} or \code{n} which 
#' has already been used in the underlying C++ code.
#' These names cannot conflict with the names of the user-defined test statistic functions 
#' \code{tsf0} and \code{tsf}.
#' The second column of \code{epList} contains the data structure 
#' of the extra parameters. 
#' Currently you can only use matirx, vector and scalor as extra parameters in 
#' the user-defined test statistic functions. 
#' A specification "matrix" is transformed to \code{arma::Mat}, 
#' and "vector" is transformed to \code{arma::Col}.
#' All other specifications are transformed to a scalor (i.e., no data structure is used).
#' The third column of \code{epList} contains the type of data in the data structure. 
#' It must be valid C++ primitive types (e.g., double, float, int and bool).
#' The following are some examples of specification of extra parameters.
#' \code{c("x", "matrix", "double")}: 1 extra parameter \code{x} with type \code{arma::Mat<double>}.
#' \code{c("x", "vector", "double")}: 1 extra parameter \code{x} with type \code{arma::Col<double>}.
#' \code{c("x", "", "int")}: 1 extra parameter \code{x} with type \code{int}.
#' 
#' @return an R function which calls C++ code for doing sequential permutation test. 
#' The parameters of the returned R function has the same meaning with parameters in 
#' the R function \code{\link{spt}}.
#' 
#' @seealso \code{\link{spt}} which is a pure R version function for doing sequential permutation test. 
#' 
#' @export
#'
#' @examples
#' # read source code
#' tsf_code = read_code("test_statistics.cpp")
#' # genrate function for sequential permutation test 
#' # for no location difference between two samples
#' if(!exists("cxxspt.mean")){
#'     cxxwrapper(tsf_code, "abs_mean_diff0", "abs_mean_diff", c("n1", "", "int")) -> cxxspt.mean
#' }
#' # download data 
#' if(!file.exists('spt_data.rda')){
#'     download.file('http://dclong.github.io/media/spt/spt_data.rda', 'spt_data.rda')
#' }
#' load('spt_data.rda')
#' cxxspt.mean(leukemia, 10, 1000, n1=5) -> cxxspt.mean.out
#' # generate function for sequential permutation test
#' # for no correlation between a response variable and a bunch of covariates
#' if(!exists("cxxspt.corr")){
#'     cxxwrapper(tsf_code, "max_abs_corr0", "max_abs_corr", c("X", "matrix", "double")) -> cxxspt.corr
#' }
#' # center and standardize columns of marker
#' X = scale(t(marker))
#' # center and standardized rows of barley 
#' Y = t(scale(t(barley)))
#' # for the purpose of illustration, only test the first 100 genes since it's time consuming
#' cxxspt.corr(Y[1:100,], 10, 1000, X=X) -> cxxspt.corr.out
#' 
cxxwrapper = function(tsf_code, tsf0, tsf, epList){
    if(missing(epList)){
        src = "arma::mat _data(Rcpp::as<arma::mat>(data));\nint _h(Rcpp::as<int>(h));\nint _n(Rcpp::as<int>(n));\n"
        src = paste(src, "\nRNGScope scope;\nreturn Rcpp::wrap(spt(_data, ", tsf0, ", ", tsf, ", _h, _n))", sep="")
        sig = 'signature(data="numeric", h="numeric", n="numeric")'
    }else{
        # parsing extra.pars
        if(is.vector(epList)){
            epList = t(epList)
        }
        apply(epList, 1, parse_parameter) -> def_args;
        src = "arma::mat _data(Rcpp::as<arma::mat>(data));\nint _h(Rcpp::as<int>(h));\nint _n(Rcpp::as<int>(n));\n"
        src = paste(src, paste(def_args[1,], collapse="\n"), sep="")
        src = paste(src, "\nRNGScope scope;\nreturn Rcpp::wrap(spt(_data, ", tsf0, ", ", tsf, ", _h, _n, ", sep="")
        src = paste(src, paste(def_args[2,], collapse=", "), "));\n", sep="")
        sig = 'signature(data="numeric", h="numeric", n="numeric", '
        sig = paste(sig, paste(paste(epList[,1], '="numeric"', sep=""), collapse=", "), ")", sep="")
    }
    includes = paste(tsf_code, "\n\n", read_code("spt.cpp"), sep="")
    settings = inline::getPlugin("RcppArmadillo")
    settings$env$PKG_CXXFLAGS = paste('-std=c++0x ',settings$env$PKG_CXXFLAGS)
    fx = inline::cxxfunction(eval(parse(text=sig)), body=src, includes=includes, 
                     plugin="RcppArmadillo",settings=settings)
    cxxspt = function(data, h, n, ...){
        p = fx(data, h, n, ...)
        r = list(p=p, h=h, n=n)
        class(r) = "spt"
        r
    }
    cxxspt
}

#' @noRd
#' Parses the parameter list into C++ code.
#' @param par a string vector with length 3. 
#' par[1]: name of the parameter
#' par[2]: data structure of the parameter 
#' par[3]: type of data
parse_parameter = function(par){
    name = par[1]
    structure = tolower(par[2])
    type = par[3]
    if(structure=="matrix"){
        type = paste(" arma::Mat<", type, "> ", sep="")
    } else if(structure=="vector"){
        type = paste(" arma::Col<", type, "> ", sep="")
    }
    def = paste(type, " _", name, "(as<", type, ">(", name, "));", sep="")
    c(def, paste("_", name, sep=""))
}

