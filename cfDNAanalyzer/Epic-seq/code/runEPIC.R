##----------------------------------------------------------------
## Copyright Stanford University 2021
## M.S.Esfahani, et al., Alizadeh and Diehn labs, Stanford University
##----------------------------------------------------------------

suppressPackageStartupMessages(library("optparse"))

# 20240711 zjp modified
# myfactors = c("bamdir","tssinfo","outdir","panelbed","targeted","groupref","groupsize","mapq","skipHist","skipOCFFrags")
# types = c("character","character","character","character","character","character","double","double","double","double")
# myfactors2 = c("b","t","o","p","a","g","r","m","s","f")
myfactors = c("epic_dir","bamdir","tssinfo","outdir","panelbed","targeted","groupref","groupsize","mapq","skipHist","skipOCFFrags")
types = c("character","character","character","character","character","character","character","double","double","double","double")
myfactors2 = c("d","b","t","o","p","a","g","r","m","s","f")
# end

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
opt = opt[which(!is.na(opt))]

runepic <- function(inputlist)
{
  # 20240711 zjp modified
  epic_dir = inputlist$epic_dir
  source(paste0(epic_dir,'/epicseq_libs.R'))
  # end 
  
  # 2024.11.24 zjp modified
  # if (is.null(inputlist$bamdir))
  # {
  #   error("A directory of bam files is required.")
  # }else{
  #   bamdir = inputlist$bamdir
  # }
  # if (is.null(inputlist$tssinfo))
  # {
  #   tssinfo = paste0("./priordata/all.tss.genes.canonical.ensembl75.txt")
  # }else{
  #   tssinfo = inputlist$tssinfo
  # }
  # if (is.null(inputlist$panelbed))
  # {
  #   panelbed = paste0(epic_dir,"/selectors_windows/all.tss.genes.canonical.ensembl75.selector.txt")
  # }else{
  #   panelbed = inputlist$panelbed
  # }
  # if (is.null(inputlist$outdir))
  # {
  #   error("An output directory must exist and provided by user.")
  # }else{
  #   outdir = inputlist$outdir
  #   dir.create(outdir, showWarnings = FALSE)
  # }
  bamdir = inputlist$bamdir
  tssinfo = inputlist$tssinfo
  panelbed = inputlist$panelbed
  outdir = inputlist$outdir
  dir.create(outdir, showWarnings = FALSE)
  # end 
  
  
  if (is.null(inputlist$mapq))
  {
    mapq = 30
  }else{
    mapq = inputlist$mapq
  }
  if (is.null(inputlist$targeted))
  {
    targeted = "wgs"
  }else{
    targeted = inputlist$targeted
    
  }
  
  if (is.null(inputlist$groupref))
  {
    groupref = "pbmc_genome"
  }else{
    groupref = inputlist$groupref
    
  }
  if (is.null(inputlist$groupsize))
  {
    groupsize = 0
  }else{
    groupsize = inputlist$groupsize
  }
  if (is.null(inputlist$skipHist))
  {
    skipHist = 0
  }else{
    skipHist = inputlist$skipHist
  }
  
  if (is.null(inputlist$skipOCFFrags))
  {
    skipOCFFrags = 0
  }else{
    skipOCFFrags = inputlist$skipOCFFrags
  }
  if (toupper(targeted)=="YES")
  {
    targeted = "epicseq"
  }else{
    targeted = "wgs"
  }
  if (groupsize<=1)
  {
    groupsize = 0
  }
  if (targeted=="epicseq"&groupsize!=0)
  {
    warning("Gene grouping is not recommended for the targeted data")
  }
  
  write("Running the shell wrapper...", stderr())
  # 20240711 zjp modified
  # epiccommand = paste("./epic_wrapper.sh ",bamdir,tssinfo,outdir,panelbed,targeted,groupsize,
  #                     groupref,mapq,skipHist,skipOCFFrags,0,0,0,sep=" ")
  epiccommand = paste(paste0(epic_dir,"/epic_wrapper.sh "),bamdir,tssinfo,outdir,panelbed,targeted,groupsize,
                      groupref,mapq,skipHist,skipOCFFrags,0,1,1,epic_dir,sep=" ")
  # end
  write("Ruuning the following command:", stderr())
  write(epiccommand, stderr())
  system(epiccommand)
}
runepic(opt)
