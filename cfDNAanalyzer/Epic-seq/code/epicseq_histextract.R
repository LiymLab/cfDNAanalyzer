##----------------------------------------------------------------
## Copyright Stanford University 2021
## M.S.Esfahani, et al., Alizadeh and Diehn labs, Stanford University
##----------------------------------------------------------------

# 20240711 zjp modified
# source('epicseq_libs.R')
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library("optparse"))
# end

if (!suppressMessages(require('optparse'))) {
  suppressMessages(install.packages('optparse',repos="http://cran.r-project.org"))
  if (!suppressMessages(require('optparse'))) {
    stop('The package optparse was not installed')
  }
}
myfactors = c("outdir","bamdir","tsspath","rnaseq","strsave","mapq","path2dep")
myfactors2 = c("o","b","t","r","s","m","p")
types = c("character","character","character","character","character","double","character")
option_list = list()
cc_cnt = 0
for (cc in myfactors)
{
  cc_cnt = cc_cnt + 1
  option_list = c(option_list, make_option(c(paste("-",myfactors2[cc_cnt],sep=""), paste("--",cc,sep="")), type = types[cc_cnt], default = NA,
                                           help = paste("This is the status of input parameters ",cc), metavar = types[cc_cnt]))
}
opt_parser = (OptionParser(option_list = option_list))
opt = (parse_args(opt_parser))
print(opt)
opt = opt[which(!is.na(opt))]

sweep_upstream <- function(inputlist)
{
  dependencies = inputlist$path2dep
  source(paste0(dependencies,'/epicseq_sidefuncs.R'))
  windows = list(up1000down1000 = c(1000,1000),
                 up1000downN750 = c(1000,-750),
                 upN750down1000 = c(-750,1000),
                 up150down50 = c(150,50))
  if (is.null(inputlist$outdir))
  {
    error("An output directory is required")
  }else{
    basepath = inputlist$outdir
  }
  
  bamdir = inputlist$bamdir
  if (is.null(inputlist$tsspath))
  {
    error("A TSS path is required")
    
  }else{
    tsspath = inputlist$tsspath
  }
  
  if (is.null(inputlist$strsave))
  {
    strsave = 'allgenes'
  }else{
    strsave = inputlist$strsave
  }
  if (is.null(inputlist$mapq))
  {
    mapq = 25
  }else{
    mapq = inputlist$mapq
  }
  
  for (tsswindow in names(windows))
  {
    uu = windows[[tsswindow]][1]
    dd = windows[[tsswindow]][2]
    mypath <- creatergnfile(tsspath,uu,dd,basepath)
    newoutdir = paste0(basepath,"/window_","upstream",uu,"_downstream_",dd)
    dir.create(newoutdir)
    selectorpath = mypath
    outputdir = newoutdir
    print(paste(paste0(dependencies,"/sizedist_parallel.sh ") ,bamdir, selectorpath, outputdir ,mapq,dependencies, sep = " "))
    system(paste(paste0(dependencies,"/sizedist_parallel.sh "),bamdir, selectorpath, outputdir , mapq, dependencies, sep = " "))
    # newoutdir2 = paste0(basepath,"/",tsswindow)
    # # dir.create(newoutdir2)
    # # system(paste0("cp ",newoutdir,"/*txt ",newoutdir2))
  }
}
sweep_upstream(opt)
