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
myfactors = c("mode","newdirocf","outdir","groupsize","tssinfopath","sortmode","skipentropy","skipocf","skipdepth","path2dep")
myfactors_short = c("m","p","o","g","t","r","u","v","w","z")
types = c("character","character","character","double","character","character","double","double","double","character")
option_list = list()
cc_cnt = 0
for (cc in myfactors)
{
  cc_cnt = cc_cnt + 1
  option_list = c(option_list, make_option(c(paste("-",substr(myfactors_short[cc_cnt],1,1),sep=""), paste("--",cc,sep="")), type = types[cc_cnt], default = NA,
                                           help = paste("This is the status of the input parameter ",cc), metavar = types[cc_cnt]))
}
opt_parser = (OptionParser(option_list = option_list))
opt = (parse_args(opt_parser))
opt = opt[which(!is.na(opt))]

generate_data <- function(inputfile)
{
  
  path2dep = inputfile$path2dep
  path2dep <<- path2dep
  source(paste0(path2dep,'/epicseq_entropymat.R'))
  source(paste0(path2dep,'/epicseq_ndrmat.R'))
  source(paste0(path2dep,'/epicseq_initialization.R'))
  source(paste0(path2dep,'/epicseq_sidefuncs.R'))
  source(paste0(path2dep,'/epicseq_ocfmat.R'))
  source(paste0(path2dep,'/epicseq_group_by_gene.R'))
  source(paste0(path2dep,'/epicseq_grouping_entities.R')) 
  mode = inputfile$mode
  newdirocf = inputfile$newdirocf
  outdir = inputfile$outdir
  newdir1 = paste0(outdir,"/up1000down1000/")
  newdir2 = paste0(outdir,"/upN750down1000/")
  newdir3 = paste0(outdir,"/up1000downN750/")
  newdir4 = paste0(outdir,"/up150down50/")
  
  groupsize = inputfile$groupsize
  tssinfopath = inputfile$tssinfopath
  sortmode = tolower(inputfile$sortmode)
  if (is.null(inputfile$skipentropy))
  {
    skipentropy = 0
  }else{
    skipentropy = inputfile$skipentropy
  }
  if (is.null(inputfile$skipocf))
  {
    skipocf = 0
  }else{
    skipocf = inputfile$skipocf
  }
  if (is.null(inputfile$skipdepth))
  {
    skipdepth = 0
  }else{
    skipdepth = inputfile$skipdepth
  }
  gepreflist = list(pbmc_exome = pbmc_exome, pbmc_genome = pbmc_genome)
  gepref = gepreflist[[sortmode]]
  tss = fread(tssinfopath,data.table=F)
  gepref = unique(gepref[gepref[,1]%in%tss[["Gene-Symbol"]],])
  gsize = groupsize
  neggens = tss[["Gene-Symbol"]][tss$Category=="negativeControl"]
  path2neg = paste0(outdir,"/negativeControl.genes.txt")
  write.table(as.matrix(neggens),path2neg,quote=F,sep="\t",row.names=F,col.names=F)
  if (groupsize>=1)
  {
    outofgroup <- groupme(gsize,gepref,gcpath,path2neg,outdir,tssinfopath)
    print("*********Wrote to DISK II**********")
    gepref = outofgroup[["gepref"]]
    gcpath = outofgroup[["gcpath"]]
    agg_gc = outofgroup[["agg_gc"]]
    overlapfile = outofgroup[["overlapfile"]]
  }
  ## generate paths:
  if (mode=="epicseq")
  {
    sigpaths1 = paste(newdir1,list.files(newdir1,pattern = "up.1000.down.1000.txt"),sep="")
    sigpaths1_bg1 = gsub("up1000down1000","upN750down1000",gsub("up.1000.down.1000",'up.-750.down.1000',sigpaths1))
    sigpaths1_bg2 = gsub("up1000down1000","up1000downN750",gsub("up.1000.down.1000",'up.1000.down.-750',sigpaths1))
    sigpaths1_ndr = gsub("up1000down1000","up150down50",gsub("up.1000.down.1000",'up.150.down.50',sigpaths1))
    sigpaths1_2k = sigpaths1
    
    names = gsub(newdir1,"", sub(".up.1000.down.1000.txt", "", sigpaths1))
    names_bg1=gsub(newdir2,"",sub(".up.-750.down.1000.txt","",sigpaths1_bg1))
    names_bg2=gsub(newdir3,"",sub(".up.1000.down.-750.txt","",sigpaths1_bg2))
    names_ndr=gsub(newdir4,"",sub(".up.150.down.50.txt","",sigpaths1_ndr))
    names_2k=gsub(newdir1,"",sub(".up.1000.down.1000.txt","",sigpaths1_2k))
    
    sigpaths1_ocf = paste(newdirocf,list.files(newdirocf,pattern = ".ocf.values.txt"),sep="")
    names_ocf = gsub(newdirocf,"",gsub(".ocf.values.txt","",sigpaths1_ocf))
    names_ocf = gsub(newdirocf,"",gsub(".ocf.values.txt","",names_ocf))
    
    sigpaths1_ocf_list = as.list(sigpaths1_ocf)
    names(sigpaths1_ocf_list) = names_ocf
    
  }else{
    sigpaths1 = paste(newdir1,list.files(newdir1,pattern = "up.1000.down.1000.txt"),sep="")
    sigpaths1_bg1 = gsub("up1000down1000","upN750down1000",gsub("up.1000.down.1000",'up.-750.down.1000',sigpaths1))
    sigpaths1_bg2 = gsub("up1000down1000","up1000downN750",gsub("up.1000.down.1000",'up.1000.down.-750',sigpaths1))
    sigpaths1_ndr = gsub("up1000down1000","up150down50",gsub("up.1000.down.1000",'up.150.down.50',sigpaths1))
    sigpaths1_2k = sigpaths1
    
    names=gsub(newdir1,"",gsub(".up.1000.down.1000.txt","",sigpaths1))
    names_bg1=gsub(newdir2,"",gsub(".up.-750.down.1000.txt","",sigpaths1_bg1))
    names_bg2=gsub(newdir3,"",gsub(".up.1000.down.-750.txt","",sigpaths1_bg2))
    names_ndr=gsub(newdir4,"",gsub(".up.150.down.50.txt","",sigpaths1_ndr))
    names_2k=gsub(newdir4,"",gsub(".up.1000.down.1000.txt","",sigpaths1_2k))
    
    names=gsub(newdir1,"",gsub(".up.1000.down.1000.txt","",names))
    names_bg1=gsub(newdir2,"",gsub(".up.-750.down.1000.txt","",names_bg1))
    names_bg2=gsub(newdir3,"",gsub(".up.1000.down.-750.txt","",names_bg2))
    names_ndr=gsub(newdir4,"",gsub(".up.150.down.50.txt","",names_ndr))
    names_2k=gsub(newdir1,"",gsub(".up.1000.down.1000.txt","",names_2k))
    
    sigpaths1_ocf = paste(newdirocf,list.files(newdirocf,pattern = ".ocf.values.txt"),sep="")
    names_ocf = gsub(newdirocf,"",gsub(".ocf.values.txt","",sigpaths1_ocf))
    sigpaths1_ocf_list = as.list(sigpaths1_ocf)
    names(sigpaths1_ocf_list) = names_ocf
    if (groupsize>1)
    {
      sigpaths1 <- group_by_gene(sigpaths1,"up.1000.down.1000.txt","hist",groupsize,outdir)
      print("grouping I is finished")
      sigpaths1_bg1 <- group_by_gene(sigpaths1_bg1,"up.-750.down.1000.txt","hist",groupsize,outdir)
      print("grouping II is finished")
      sigpaths1_bg2 <- group_by_gene(sigpaths1_bg2,"up.1000.down.-750.txt","hist",groupsize,outdir)
      print("grouping III is finished")
      sigpaths1_ndr <- group_by_gene(sigpaths1_ndr,"up.150.down.50.txt","hist",groupsize,outdir)
      print("grouping IV is finished")
      sigpaths1_2k = sigpaths1
      sigpaths1_ocf <- group_by_gene(sigpaths1_ocf,".ocf.values.txt","ocf",groupsize,outdir)
      print("grouping of OCF is finished")
      sigpaths1_ocf_list = as.list(sigpaths1_ocf)
      names(sigpaths1_ocf_list) = names_ocf
      print("*********Gene Grouping is finished**********")
    }
  }
  sigpaths1_list = as.list(sigpaths1)
  names(sigpaths1_list) = names
  sigpaths1_bg1_list = as.list(sigpaths1_bg1)
  names(sigpaths1_bg1_list) = names_bg1
  sigpaths1_bg2_list = as.list(sigpaths1_bg2)
  names(sigpaths1_bg2_list) = names_bg2
  ndr_list = as.list(sigpaths1_ndr)
  names(ndr_list) = names_ndr
  ndr2k_list = as.list(sigpaths1_2k)
  names(ndr2k_list) = names_2k
  bgs_list = list()
  for (nm in names)
  {
    bgs_list[[nm]] = c(sigpaths1_bg1_list[[nm]],sigpaths1_bg2_list[[nm]])
  }
  allsigs = sigpaths1_list
  allndr = ndr_list
  all2k = ndr2k_list
  allbgs = bgs_list
  allocf = sigpaths1_ocf_list
  #### Entropy
  if (skipentropy==0)
  {  
    entropy_matrix = epic_calcpfe(allsigs,allbgs, gsize=groupsize,mode)
    entropyfile_path = paste(outdir,"/pfe.matrix.merged.by.",groupsize,".txt",sep="")
    write.table(entropy_matrix,entropyfile_path,quote=F,sep="\t",row.names=F)
  }
  print("*********PFE is calculated**********")
  
  #### NDR Coverage:
  if (skipdepth==0)
  { 
    ndr_matrix = epic_calcndr(allndr, gsize=groupsize)
    ndrfile_path = paste(outdir,"/coverage.ndr.matrix.merged.by.",groupsize,".txt",sep="")
    write.table(ndr_matrix,ndrfile_path,quote=F,sep="\t",row.names=F)
    #### 2k Coverage:
    d2k_matrix = epic_calcndr(all2k, gsize=groupsize)
    d2kfile_path = paste(outdir,"/coverage.2k.matrix.merged.by.",groupsize,".txt",sep="")
    write.table(d2k_matrix,d2kfile_path,quote=F,sep="\t",row.names=F)
    #### Flanking Coverage:
    flanking_matrix = epic_calcndr(allbgs, gsize=groupsize)
    flankingfile_path = paste(outdir,"/coverage.flanking.matrix.merged.by.",groupsize,".txt",sep="")
    write.table(flanking_matrix,flankingfile_path,quote=F,sep="\t",row.names=F)
  }
  #### OCF 
  if (skipocf==0)
  { 
    print("****************************")
    print("working on OCF Data")
    ocf_matrix = epic_calcocf(allocf, gsize=groupsize)
    ocffile_path = paste(outdir,"/ocf.matrix.merged.by.",groupsize,".txt",sep="")
    write.table(ocf_matrix,ocffile_path,quote=F,sep="\t",row.names=F)
  } 
}
generate_data(opt)
