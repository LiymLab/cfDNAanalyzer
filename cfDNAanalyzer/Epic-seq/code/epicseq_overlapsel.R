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

myfactors = c("tsspath","selector","outdir","path2dep")
myfactors2 = c("t","s","o","p")
types = c("character","character","character","character")
option_list = list()
cc_cnt = 0
for (cc in myfactors)
{
  cc_cnt = cc_cnt + 1
  option_list = c(option_list, make_option(c(paste("-",myfactors2[cc_cnt],sep=""), paste("--",cc,sep="")), type = types[cc_cnt], default = NA,
                                           help = paste("This is the status of input parameter ",cc), metavar = types[cc_cnt]))
}
opt_parser = (OptionParser(option_list = option_list))
opt = (parse_args(opt_parser))
opt = opt[which(!is.na(opt))]
calc_overlap <- function(inputlist)
{
  path2dep = inputlist$path2dep
  print(path2dep)
  source(paste0(path2dep,"/epicseq_sidefuncs.R"))
  tsspathpath = inputlist$tsspath
  selector = inputlist$selector
  outdir = inputlist$outdir
  #### re-creating the selector files
  path1 = creatergnfile(tsspathpath, uu = 150, dd = 50,outdir)
  path2 = creatergnfile(tsspathpath, uu = 1000, dd = 1000,outdir)
  path3 = creatergnfile(tsspathpath,uu=1000,dd=-750,outdir)
  path4 = creatergnfile(tsspathpath,uu=-750,dd=1000,outdir)
  
  ll1 = fread(path1,data.table=F)
  ll2 = fread(path2,data.table=F)
  tsscords = fread(tsspathpath,data.table=F)
  ll1$lookup = paste(tsscords[,1],tsscords[,2],tsscords[,3],sep="@@")
  ll2$lookup = paste(tsscords[,1],tsscords[,2],tsscords[,3],sep="@@")
  
  write.table(ll1,paste(outdir,"/temp.ndr.150.50.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
  write.table(ll2,paste(outdir,"/temp.ndr.1000.1000.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
  path1new = paste(outdir,"/temp.ndr.150.50.txt",sep="")
  path2new = paste(outdir,"/temp.ndr.1000.1000.txt",sep="")
  
  overlapfile1 = paste(outdir,"/overlap.selector.150.50.txt",sep="")
  overlapfile2 = paste(outdir,"/overlap.selector.1000.1000.txt",sep="")
  intersectstr1 = paste("bedtools intersect -a ",path1new," -b ",selector," -wo | sort -u >", overlapfile1,sep="")
  intersectstr2 = paste("bedtools intersect -a ",path2new," -b ",selector," -wo | sort -u >", overlapfile2,sep="")
  system(intersectstr1)
  system(intersectstr2)
  ovl_data1 = fread(overlapfile1,data.table=F,header = F)
  ovl_data2 = fread(overlapfile2,data.table=F,header = F)
  ovl_data1$lookupToCollapse1 = paste(ovl_data1[,1],ovl_data1[,2],ovl_data1[,3],ovl_data1[,4],sep="@@")
  ovl_data1$lookup = ovl_data1[,7]
  ovl_data2$lookupToCollapse2 = paste(ovl_data2[,1],ovl_data2[,2],ovl_data2[,3],ovl_data2[,4],sep="@@")
  ovl_data2$lookup = ovl_data2[,7]
  
  ovl_data1_uniq = unique(ovl_data1[c("lookupToCollapse1","lookup")])
  ovl_data2_uniq = unique(ovl_data2[c("lookupToCollapse2","lookup")])
  
  aggregated1 = aggregate(x=ovl_data1[c("V10")],by=list(lookup = ovl_data1[["lookup"]]),FUN=sum)
  aggregated1$V10 = pmin(aggregated1$V10,200)
  aggregated2 = aggregate(x=ovl_data2[c("V10")],by=list(lookup = ovl_data2[["lookup"]]),FUN=sum)
  aggregated2$V10 = pmin(aggregated2$V10,2000)
  
  colnames(aggregated1) = c("lookup","SelectorNDRbases")
  aggregated1 = merge(aggregated1,ovl_data1_uniq,by="lookup")
  colnames(aggregated2) = c("lookup","Selector2Kbases")
  aggregated2 = merge(aggregated2,ovl_data2_uniq,by="lookup")
  aggregate_df = merge(aggregated1,aggregated2,by="lookup")
  pp = ((unlist(apply(aggregate_df,1,function(x){y = x[["lookupToCollapse2"]];z=as.character(y);u=unlist(strsplit(z,"@@"));return(u)}))))
  pp = t(matrix(pp,nrow=4))
  pp = as.data.frame(pp)
  colnames(pp) = c("chr","start","end","genename")
  myfull_df = cbind(pp,aggregate_df[c("SelectorNDRbases","Selector2Kbases")])
  myfull_df = renameTSS(myfull_df)
  tsspathname = basename(tsspathpath)
  nn = gsub(".txt","",tsspathname)
  file2write = paste(outdir,"/",nn,".overlaping.bases.info.file.txt",sep="")
  write.table(myfull_df,file2write,quote=F,sep="\t",row.names=F)
  return(file2write)
}
calc_overlap(opt)

