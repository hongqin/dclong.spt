#!/usr/bin/Rscript
library(roxygen2)
roxygenize('dclong.spt')
f = 'dclong.spt/R/.Rhistory'
if(file.exists(f)){
    file.remove(f)
}

