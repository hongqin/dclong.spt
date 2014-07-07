#' @title Read in C++ code 
#' @description The function \code{read_code} read in C++ code 
#' in the "rcpp" folder under the home directory of this package.
#' @param title title of the file containing (C++) code.
#' @examples
#' read_code('spt.cpp')
#' read_code('test_statistics.cpp')
#' @export
read_code = function(title){
    system.file(paste("rcpp/",title,sep=""), package="dclong.spt") -> path
    paste(readLines(con=path), collapse="\n")
}
