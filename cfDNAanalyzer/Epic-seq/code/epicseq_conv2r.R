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


myfactors = c("tssinfo","outdir","mode","path2dep")
myfactors2 = c("t","o","m","p")
types = c("character","character","character","character","character")
option_list = list()
cc_cnt = 0
for (cc in myfactors)
{
  cc_cnt = cc_cnt + 1
  option_list = c(option_list, make_option(c(paste("--",cc,sep="")), type = types[cc_cnt], default = NA,
                                           help = paste("This is the status of optimization factor ",cc), metavar = types[cc_cnt]))
}
opt_parser = (OptionParser(option_list = option_list))
opt = (parse_args(opt_parser))
opt = opt[which(!is.na(opt))]
convert2readable <- function(inputfile)
{
  write("Converting samtools-based histograms into matrices...", stderr())
  path2dep = inputfile$path2dep
  source(paste0(path2dep,"/epicseq_featuremat.R"))
  source(paste0(path2dep,"/epicseq_featuremat_ocf.R"))
  tssinfo = inputfile$tssinfo
  outdir = inputfile$outdir
  mode = inputfile$mode
  ocfout = paste(outdir,"/ocfoutputs",sep="")
  subdirs = c("window_upstream1000_downstream_1000",
              "window_upstream1000_downstream_-750",
              "window_upstream-750_downstream_1000",
              "window_upstream150_downstream_50")
  windows = list("window_upstream1000_downstream_1000"="up1000down1000",
                 "window_upstream1000_downstream_-750" = "up1000downN750",
                 "window_upstream-750_downstream_1000" = "upN750down1000",
                 "window_upstream150_downstream_50"= "up150down50")
  p1 = gsub(".txt","_upstream_1000_downstream_1000.txt",tssinfo)
  p1 = paste(outdir,"/",basename(p1),sep="")
  p2 = gsub(".txt","_upstream_1000_downstream_-750.txt",tssinfo)
  p2 = paste(outdir,"/",basename(p2),sep="")
  p3 = gsub(".txt","_upstream_-750_downstream_1000.txt",tssinfo)
  p3 = paste(outdir,"/",basename(p3),sep="")
  p4 = gsub(".txt","_upstream_150_downstream_50.txt",tssinfo)
  p4 = paste(outdir,"/",basename(p4),sep="")
  
  selectors = c(p1,p2,p3,p4)
  
  names = c("up.1000.down.1000","up.1000.down.-750","up.-750.down.1000","up.150.down.50")
  myfullgenepath = tssinfo
  currpath = outdir
  ddcnt = 0
  for (dd in subdirs)
  {
    ddcnt = ddcnt + 1
    saveflag_str = names[ddcnt]
    path2selector = selectors[ddcnt]
    finald = paste0(currpath,"/",dd)
    create_ftr_matrix(pathtodata= finald,path2selector= path2selector,
                           saveflag_str=saveflag_str,mode=mode,myfullgenepath= myfullgenepath)
    newoutdir2 = paste0(currpath,"/",windows[[dd]])
    dir.create(newoutdir2)
    system(paste0("cp ",finald,"/*txt ",newoutdir2))
  }
  create_ftr_matrix_ocf(pathtodata = ocfout,saveflag_str = 'ocf.values',mode = mode,myfullgenepath=tssinfo)
}
convert2readable(opt)
