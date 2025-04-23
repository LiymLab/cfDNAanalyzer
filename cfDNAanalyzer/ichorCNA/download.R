# git clone git@github.com:broadinstitute/ichorCNA.git  
## install from CRAN
install.packages("plyr") 
## install packages from
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("HMMcopy")  
library("HMMcopy")
BiocManager::install("GenomeInfoDb")  
library("GenomeInfoDb")
BiocManager::install("GenomicRanges")  
library("GenomicRanges")
## from the command line and in the directory where ichorCNA github was cloned.
#R CMD INSTALL ichorCNA    

# HMMcopy (for readCounter)
# https://github.com/shahcompbio/hmmcopy_utils