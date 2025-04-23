##----------------------------------------------------------------
## Copyright Stanford University 2021
## M.S.Esfahani, et al., Alizadeh and Diehn labs, Stanford University
##----------------------------------------------------------------
source('epicseq_sidefuncs.R')
source('epicseq_libs.R')
myfactors = c("xnewpath","tssinfo","outdir")
myfactors2 = c("x","t","o")
types = c("character","character","character")
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

infer_gep <- function(inputlist, mode = "wgs")
{
  
  path2models = paste0("pretrainedmodel/",list.files("pretrainedmodel/",pattern=".rds"))
  
  if (is.null(inputlist$xnew))
  {
    error("An xnewpath should be provided.")
  }else{
    xnewpath = inputlist$xnewpath
  }
  if (is.null(inputlist$outdir))
  {
    error("An output directory should be provided.")
  }else{
    outdir = inputlist$outdir
  }
  if (is.null(inputlist$tssinfo))
  {
    error("A TSS information file (corresponding to TSSs targeted) should be provided. This file will be used for normalizations.")
  }else{
    tssinfo = inputlist$tssinfo
  }
  tssinfodf = fread(tssinfo,data.table=F)
  colnames(tssinfodf) = c("CHR","TSS","Gene-Symbol","Category","Strand")
  neggenes = tssinfodf[["Gene-Symbol"]][tssinfodf$Category=="negativeControl"]
  xnew = fread(xnewpath,data.table=F)
  ftrs = c("PFE","NDR_relative_2k")
  xnewnormalized = c()
  allneggenesmeans = c()
  for (kk in unique(xnew$Sample))
  {
    currmat = xnew[xnew$Sample==kk,]
    currmat = samplewise_adjust_final(currmat,neggenes)
    xnewnormalized = rbind(xnewnormalized,currmat)
  }
  xnewnormalized_intercepted = cbind(data.frame(intercept=rep(1,nrow(xnewnormalized))),xnewnormalized[ftrs])
  models = readRDS(path2models)
  models = models[models$mode=="combo",]
  models_coeffs = models[c("intercept",ftrs)]
  all_scores = apply(xnewnormalized_intercepted[c("intercept",ftrs)],1,
                     function(x){xreps = rep.row(x,nrow(models_coeffs));median(rowSums(xreps*models_coeffs))})
  fullscores_df = data.frame(all_scores)
  colnames(fullscores_df) = "inferredGEP"
  finaldf = cbind(xnewnormalized,fullscores_df)
  write.table(finaldf,paste0(outdir,"/EPICSeqInferredExpressionValues.txt"),quote=F,sep="\t",row.names = F)
}
infer_gep(opt)