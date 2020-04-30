rm(list=ls()) #clear out environment
library(GEOquery)  #require package


gset <- getGEO("GSE42414", GSEMatrix=TRUE) #load the NCBI database data using internet
