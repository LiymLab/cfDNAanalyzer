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

myfactors = c("outdir","tssinfo","groupsize")
types = c("character","character","double")
option_list = list()
cc_cnt = 0
for (cc in myfactors)
{
  cc_cnt = cc_cnt + 1
  option_list = c(option_list, make_option(c(paste("-",substr(cc,1,1),sep=""), paste("--",cc,sep="")), type = types[cc_cnt], default = NA,
                                           help = paste("This is the status of input parameter ",cc), metavar = types[cc_cnt]))
}
opt_parser = (OptionParser(option_list = option_list))
opt = (parse_args(opt_parser))
opt = opt[which(!is.na(opt))]
featurize_data <- function(inputfile)
{
  outdir = inputfile$outdir
  tssinfo = inputfile$tssinfo
  gsize = inputfile$groupsize
  entpath = paste(outdir,"/pfe.matrix.merged.by.",gsize,".txt",sep="")
  cov2kpath = paste(outdir,"/coverage.2k.matrix.merged.by.",gsize,".txt",sep="")
  covNDRpath = paste(outdir,"/coverage.ndr.matrix.merged.by.",gsize,".txt",sep="")
  covflankingpath = paste(outdir,"/coverage.flanking.matrix.merged.by.",gsize,".txt",sep="") 
  ocfpath = paste(outdir,"/ocf.matrix.merged.by.",gsize,".txt",sep="")
  tssinfoname = basename(tssinfo)
  nn = gsub(".txt","",tssinfoname)
  file2write = paste(outdir,"/",nn,".overlaping.bases.info.file.txt",sep="")
  if (gsize>1)
  {
    file2write = gsub(".txt",paste(".merged.by.",gsize,".txt",sep=""),file2write)
  }
  ratiocovered = fread(file2write,data.table=F)
  ########## Reading PFE File
  ## Here we first read the entropy file, and only focus on the "_Entropy" features
  print("Working on PFE")
  entropy_data = fread(entpath,data.table=F)
  print(colnames(entropy_data))
  entropy_data_limited = entropy_data[,c(1,2,3,grep("_PFE$",colnames(entropy_data)))]
  entropy_avg_data_limited = entropy_data[,c(1,2,3,grep("_Shannon$",colnames(entropy_data)))]
  entropy_gini = entropy_data[,c(1,2,3,grep("_Gini$",colnames(entropy_data)))]
  ########## Reading OCF Files & the corresponding calculations
  ## Here we calculate two new features-- Up (OCF pvalue) and Down (OCF pvalue)
  print("Working on the OCF Files")
  ocf_data = fread(ocfpath,data.table=F)
  ocf_data_RDiff = ocf_data[,c(1,2,3,grep("OCF_RDiff$",colnames(ocf_data)))]
  ocf_data_LDiff = ocf_data[,c(1,2,3,grep("OCF_LDiff$",colnames(ocf_data)))]
  ocf_data_RTot = ocf_data[,c(1,2,3,grep("OCF_RTot$",colnames(ocf_data)))]
  ocf_data_LTot = ocf_data[,c(1,2,3,grep("OCF_LTot$",colnames(ocf_data)))]
  ocf_data_signal = ocf_data_RDiff
  ocf_data_signal[,-c(1,2,3)] = 
    (ocf_data_RDiff[,-c(1,2,3)]-ocf_data_LDiff[,-c(1,2,3)])/
    (ocf_data_LTot[,-c(1,2,3)]+ocf_data_RTot[,-c(1,2,3)])
  colnames(ocf_data_signal)[-c(1,2,3)] = 
    gsub("OCF_RDiff","OCF_Ratio",colnames(ocf_data_signal)[-c(1,2,3)])
  ocf_data_signal_raw = ocf_data_RDiff 
  ocf_data_signal_raw[,-c(1,2,3)] = 
    ocf_data_RDiff[,-c(1,2,3)]-ocf_data_LDiff[,-c(1,2,3)]
  colnames(ocf_data_signal_raw)[-c(1,2,3)] = 
    gsub("OCF_RDiff","OCF_Raw",colnames(ocf_data_signal_raw)[-c(1,2,3)])
  ########## Reading Coverage files & the corresponding calculations
  print("Working on Coverage Files")
  cov2k_data = fread(cov2kpath,data.table=F)
  covndr_data = fread(covNDRpath,data.table=F)
  covflanking_data = fread(covflankingpath,data.table=F)  
  cov2k_data_middle = cov2k_data[,c(1,2,3,grep("NDR_middle$",colnames(cov2k_data)))]
  colnames(cov2k_data_middle)[grep("NDR_middle$",colnames(cov2k_data_middle))] = 
    gsub("NDR","D2K",colnames(cov2k_data_middle)[grep("NDR_middle$",colnames(cov2k_data_middle))])
  covNDR_data_middle = covndr_data[,c(1,2,3,grep("NDR_middle$",colnames(covndr_data)))]
  covNDR_data_middle_pval = covNDR_data_middle
  covNDR_data_middle_pois = covNDR_data_middle
  covNDR_data_middle_cpm = covNDR_data_middle
  cov2k_data_middle_cpm = cov2k_data_middle
  covNDR_data_middle_relative = covNDR_data_middle
  n = ncol(covNDR_data_middle) - 3
  nr = nrow(covNDR_data_middle) # Number of TSS regions
  tss_ids= cov2k_data_middle[,1]
  for (i in seq(n))
  {
    xndr = covNDR_data_middle[,i+3]
    x2k = cov2k_data_middle[,i+3]
    xndr_CPM = xndr/sum(xndr)*1E6
    x2k_CPM = x2k/sum(x2k)*1E6
    cov2k_data_middle_cpm[,i+3] = x2k_CPM
    covNDR_data_middle_cpm[,i+3] = xndr_CPM
  }
  colnames(covNDR_data_middle_cpm)[-c(1,2,3)] = 
    gsub("_middle","_middle_cpm",colnames(covNDR_data_middle_cpm)[-c(1,2,3)])
  cov2k_data_all = cov2k_data[,c(1,2,3,grep("NDR_all$",colnames(cov2k_data)))]
  colnames(cov2k_data_all)[grep("NDR_all$",colnames(cov2k_data_all))] = 
    gsub("NDR","D2K",colnames(cov2k_data_all)[grep("NDR_all$",colnames(cov2k_data_all))])
  covflanking_data_all = covflanking_data[,c(1,2,3,grep("NDR_all$",colnames(covflanking_data)))]
  colnames(covflanking_data_all)[grep("NDR_all$",colnames(covflanking_data_all))] =
    gsub("NDR","FLNK",colnames(covflanking_data_all)[grep("NDR_all$",colnames(covflanking_data_all))])
  covNDR_data_all = covndr_data[,c(1,2,3,grep("NDR_all$",colnames(covndr_data)))]
  covNDR_data_all_pval = covNDR_data_all
  covNDR_data_all_pois = covNDR_data_all
  covNDR_data_all_cpm = covNDR_data_all
  cov2k_data_all_cpm = cov2k_data_all
  covNDR_data_all_relative = covNDR_data_all
  covflanking_data_all_cpm = covflanking_data_all
  n = ncol(covNDR_data_all) - 3
  nr = nrow(covNDR_data_all) 
  tss_ids= cov2k_data_all[,1]
  for (i in seq(n))
  {
    xndr = covNDR_data_all[,i+3]
    x2k = cov2k_data_all[,i+3]
    xflnk = covflanking_data_all[,i+3]
    xndr_CPM = xndr/sum(xndr)*1E6
    x2k_CPM = x2k/sum(x2k)*1E6
    xflnk_CPM = xflnk/sum(xflnk)*1E6
    cov2k_data_all_cpm[,i+3] = x2k_CPM
    covNDR_data_all_cpm[,i+3] = xndr_CPM
    covflanking_data_all_cpm[,i+3]= xflnk_CPM
  }
  colnames(covNDR_data_all_cpm)[-c(1,2,3)] = 
    gsub("_all","_all_cpm",colnames(covNDR_data_all_cpm)[-c(1,2,3)])
  ########## Merging all features across all TSSs and samples in the cohort---
  ## Next we combine all the features (5-10??) for all the TSS's and samples--
  ## This is very similar to a MAF file- where each row corresponds to one sample and one TSS
  write("Working on combining features to create the meta file", stderr())
  final_ftrs = merge(covNDR_data_middle_cpm,cov2k_data_middle_cpm,
                     by=c("TSS_ID","gc","gene") )
  final_ftrs = merge(final_ftrs,cov2k_data_all_cpm,by=c("TSS_ID","gc","gene") )
  final_ftrs = merge(final_ftrs,covNDR_data_all_cpm,by=c("TSS_ID","gc","gene") )
  final_ftrs = merge(final_ftrs,covflanking_data_all_cpm,by=c("TSS_ID","gc","gene") )
  final_ftrs = merge(final_ftrs,entropy_data_limited,by=c("TSS_ID","gc","gene") )
  final_ftrs = merge(final_ftrs,entropy_avg_data_limited,by=c("TSS_ID","gc","gene") )
  final_ftrs = merge(final_ftrs,entropy_gini,by=c("TSS_ID","gc","gene") )
  final_ftrs = merge(final_ftrs,ocf_data_signal,by=c("TSS_ID","gc","gene") )
  final_ftrs = merge(final_ftrs,ocf_data_signal_raw,by=c("TSS_ID","gc","gene") )
  
  featurestrs = c("PFE","Shannon","Gini",
                  "D2K_all","NDR_all_cpm","OCF_Ratio")
  allftrsmelted = c()
  ftr_cnt = 0
  for (ftr in featurestrs)
  {
    ftr_cnt = ftr_cnt + 1
    currftrs = final_ftrs[,c(1,2,3,grep(ftr,colnames(final_ftrs)))]
    colnames(currftrs) = gsub(paste("_",ftr,sep=""),"",colnames(currftrs))
    ui=melt(currftrs,id.vars=c("TSS_ID","gc","gene"))
    colnames(ui)[4] = "Sample"
    colnames(ui)[5] = ftr
    if (ftr_cnt==1)
    {
      allftrsmelted = ui
    }else{
      allftrsmelted = merge(allftrsmelted,ui,by=c("TSS_ID","gc","gene","Sample"))
    }
  }
  allftrsmelted$NDR_relative_2k = allftrsmelted$NDR_all_cpm/(0.001+allftrsmelted$D2K_all)
  allftrsmelted$NDR_all_cpm = NULL
  allftrsmelted$D2K_all = NULL
  write.table(allftrsmelted,paste(outdir,'/epicseq.features.merged.by.',gsize,'.txt',sep=""),quote=F,sep="\t",row.names=F)
  write("***** All files are created *****", stderr())
  
}
featurize_data(opt)