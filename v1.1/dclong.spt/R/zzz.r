#' @noRd
#' Load the dynamic library when this library is loaded.
.onLoad <- function(lib,pkg){
  library.dynam("dclong.spt",pkg,lib)  
}
