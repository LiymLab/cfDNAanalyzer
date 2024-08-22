##----------------------------------------------------------------
## Copyright Stanford University 2021
## M.S.Esfahani, et al., Alizadeh and Diehn labs, Stanford University
##----------------------------------------------------------------

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library("optparse"))
if (!suppressMessages(require('optparse'))) {
  suppressMessages(install.packages('optparse',repos="http://cran.r-project.org"))
  if (!suppressMessages(require('optparse'))) {
    stop('The package optparse was not installed')
  }
}