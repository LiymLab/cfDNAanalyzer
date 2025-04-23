##----------------------------------------------------------------
## Copyright Stanford University 2021
## M.S.Esfahani, et al., Alizadeh and Diehn labs, Stanford University
##----------------------------------------------------------------
create_ftr_matrix <- function(pathtodata, path2selector,saveflag_str= NULL, mode = "wgs", myfullgenepath=NULL)
{
  mode = tolower(mode)
  if (is.null(myfullgenepath)){
    error("if mode is not WGS or EPICSEQ then user must enter a path to the full info of the file")
  }else{
    myfullgene = fread(myfullgenepath,data.table=F)
  }
  colnames(myfullgene) = c("CHR","TSS","Gene-Symbol","Category","Strand")
  myfullgene_unique = unique(myfullgene[c("CHR","TSS","Gene-Symbol","Category","Strand")])
  basepath = pathtodata
  annotation_file = fread(path2selector,data.table = F)
  annotation_file = unique(annotation_file)
  annotation_file$V7 = myfullgene_unique$Strand[match(annotation_file$V4,myfullgene_unique[["Gene-Symbol"]])]
  colnames(annotation_file) = c("chr","start","end","genename","reason","tss","strand")
  mydirs = list.dirs(path=basepath,recursive=FALSE)  
  bins = seq(50, 398, 1)
  bin_names = paste("bin_",bins[seq(1,length(bins))], sep="")
  all_samples_gene = c()
  dr_cnt = 0
  for (dr in mydirs)
  {
    dr_cnt = dr_cnt + 1
    DR = dr
    file_bins = list.files(path = DR, pattern = "bin_out")
    for (file_bin in file_bins)
    {
      curr_info = fread(paste(DR, "/",file_bin[grep("bin_out",file_bin)], sep=""), data.table=F,sep=" ")
      curr_info = cbind( t(matrix(unlist(lapply(curr_info$V1, function(x){z= as.character(x);
      yy = unlist(strsplit(z, "\t"))})),nrow = 5)), curr_info[,-c(1)])
      curr_info = data.frame(curr_info)
      colnames(curr_info) = c("chr","start","end","gene",bin_names)
      if (length(grep("chr",curr_info$chr))==0){
        curr_info$chr = paste("chr",curr_info$chr,sep="")
      }
      curr_info_names = paste(curr_info$chr,curr_info$start,curr_info$end,curr_info$gene,sep="_")
      annotation_file_names = paste(annotation_file$chr,annotation_file$start,annotation_file$end,annotation_file$genename,sep="_")
      myoverlap = annotation_file_names
      annotation_file = annotation_file[match(myoverlap,annotation_file_names),]
      curr_info = curr_info[match(myoverlap,curr_info_names),]
      curr_info_merged = cbind(curr_info,annotation_file[c("genename","reason","strand","tss")] )
      if (sum(is.na(curr_info_merged[,1]))>0){
        curr_info_merged[is.na(curr_info_merged$chr),c(1,2,3,4)] =  matrix(annotation_file[is.na(curr_info_merged$chr),c("chr","tss","tss","genename")],ncol=4)
      }
      curr_info_merged$sampleName = gsub(".dualindex-deduped.sorted.txt","",
                                         gsub("bin_out_Sample_","",gsub(".singleindex-deduped.sorted.txt","",file_bin[grep("bin_out",file_bin)])))
      curr_info_merged$sampleName = gsub("_cfDNA.sorted.samtools-deduped.sorted.txt","",curr_info_merged$sampleName)
      curr_info_merged_sub = curr_info_merged
      if (!is.null(saveflag_str))
      {
        nn = dr
        write.table(curr_info_merged_sub,paste(nn,".",saveflag_str,".txt",sep=""),quote=F,sep="\t",row.names=F)
        
      }
    }
  }
  return()
}
