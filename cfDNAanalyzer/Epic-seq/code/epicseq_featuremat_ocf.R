##----------------------------------------------------------------
## Copyright Stanford University 2021
## M.S.Esfahani, et al., Alizadeh and Diehn labs, Stanford University
##----------------------------------------------------------------

create_ftr_matrix_ocf <- function(pathtodata,saveflag_str= NULL, myfullgenepath= NULL, mode = "wgs")
{
  myfullgene = fread(myfullgenepath,data.table=F)
  basepath = pathtodata
  colnames(myfullgene) = c("CHR","TSS","Gene-Symbol","Category","Strand")
  myfullgene_unique = unique(myfullgene[c("CHR","TSS","Gene-Symbol","Category","Strand")])
  myfullgene_unique$start = myfullgene_unique$TSS - 1000
  myfullgene_unique$end = myfullgene_unique$TSS + 1000
  myfullgene_unique$chr = myfullgene_unique$CHR
  myfullgene_unique$genename = myfullgene_unique[["Gene-Symbol"]]
  myfullgene_unique$reason = myfullgene_unique$Category
  myfullgene_unique$strand = myfullgene_unique$Strand
  annotation_file = myfullgene_unique
  annotation_file = annotation_file[c("chr","start","end","genename","reason","strand","TSS")]
  mydirs=list.dirs(path=basepath,recursive=FALSE)
  colNames = c("OCF_1","OCF_x","OCF_LDiff","OCF_LTot","OCF_RDiff","OCF_RTot")
  all_samples_gene = c()
  dr_cnt = 0
  for (dr in mydirs)
  {
    dr_cnt = dr_cnt + 1
    DR = dr
    file_bins = list.files(path = DR, pattern = "bin_out")
    for (file_bin in file_bins)
    {
      curr_info = fread(paste(DR, "/",file_bin[grep("bin_out",file_bin)], sep=""), data.table=F, sep=" ")
      curr_info = cbind( t(matrix(unlist(lapply(curr_info$V1, function(x){z= as.character(x);
      yy = unlist(strsplit(z, "\t"))})),nrow = 5)), curr_info[,-c(1)])
      curr_info = data.frame(curr_info)
      colnames(curr_info) = c("chr","TSS","genename","strand",colNames)
      if (length(grep("chr",curr_info$chr))==0){
        curr_info$chr = paste("chr",curr_info$chr,sep="")
      }
      
      curr_info_names = paste(curr_info$chr,curr_info$TSS,curr_info$strand,sep="_")
      annotation_file_names = paste(annotation_file$chr,annotation_file$TSS,annotation_file$strand,sep="_")
      myoverlap = annotation_file_names
      annotation_file = annotation_file[match(myoverlap,annotation_file_names),]
      naones = which(is.na(match(myoverlap,curr_info_names)))
      curr_info = curr_info[match(myoverlap,curr_info_names),]
      curr_info_merged = cbind(curr_info,annotation_file[c("genename","reason","strand")] )
      curr_info_merged[,1] = annotation_file$chr
      curr_info_merged[,2] = annotation_file$TSS
      curr_info_merged[,3] = annotation_file$genename
      curr_info_merged[,4] = annotation_file$reason
      curr_info_merged$sampleName = gsub(".dualindex-deduped.sorted.txt","",
                                         gsub("bin_out_Sample_","",gsub(".singleindex-deduped.sorted.txt","",file_bin[grep("bin_out",file_bin)])))
      curr_info_merged$sampleName = gsub("_cfDNA.sorted.samtools-deduped.sorted.txt","",curr_info_merged$sampleName)
      curr_info_merged_sub = curr_info_merged
      curr_ids = paste(curr_info_merged_sub[,1],curr_info_merged_sub[,2],curr_info_merged_sub[,3],sep="_")
      if (!is.null(saveflag_str))
      {
        nn = dr
        write.table(curr_info_merged_sub,paste(nn,".",saveflag_str,".txt",sep=""),quote=F,sep="\t",row.names=F)
        
      }
      old_ids = paste(all_samples_gene[,1],all_samples_gene[,2],all_samples_gene[,3],sep="_")
    }
  }
  return()
}
