#!/usr/local/bin Rscript
# test if there is at least one argument: if not, return an error
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
    print("Using default output file type")
    library(rmarkdown)
    rmarkdown::render(args[1], "pdf_document")
} else if (length(args)==2) {
    print("Using 2nd argument as output file type")
    library(rmarkdown)
    rmarkdown::render(args[1], args[2])
}
