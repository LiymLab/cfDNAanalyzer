##----------------------------------------------------------------
## Copyright Stanford University 2021
## M.S.Esfahani, et al., Alizadeh and Diehn labs, Stanford University
##----------------------------------------------------------------

groupme <- function(gsize,gepref,gcpath,path2neg,outdir,tssinfopath)
{
  gcpathfilename = basename(gcpath)
  gcpathnew = paste0(outdir,"/",gsub(".bed",paste(".grouped.by.",gsize,".bed",sep=""),gcpathfilename))
  negconts = read.table(path2neg,header=F,sep="\t")
  neggenes = as.matrix(negconts[,1])
  gcinfo = fread(gcpath,data.table=F)
  gepref_here = gepref
  gepref_here = gepref_here[gepref_here$geneName%in%gcinfo[,1],]
  gcinfo = gcinfo[gcinfo[,1]%in%gepref_here$geneName,]
  mypaste = function(x){return(paste(x,collapse="_&&_"))}
  gepref_negcont = gepref_here[gepref_here$geneName%in%neggenes,]
  gepref_rest = gepref_here[!gepref_here$geneName%in%neggenes,]
  gepref_sorted = gepref_rest[order(-gepref_rest$avg_exp),]
  avgexp = rollapply(gepref_sorted$avg_exp,width=gsize,by=gsize,FUN = mean)
  genes = rollapply(gepref_sorted$geneName,width=gsize,by=gsize,FUN = mypaste)
  neggenes_data = gepref_negcont[c("geneName","avg_exp")]
  other_genes_data = data.frame(geneName = genes, avg_exp = avgexp)
  grouped_info = rbind(other_genes_data,neggenes_data)
  write.table(grouped_info,paste(outdir,"/response.expression.grouped.by.",gsize,".txt",sep=""),quote=F,sep="\t",row.names=F)
  #############
  groups = as.matrix(grouped_info$geneName)
  gcsnew = c()
  for (g in groups)
  {
    geneshere = unlist(strsplit(as.character(g),"_&&_"))
    gccurr = gcinfo[match(geneshere,gcinfo[,1]),]
    gcsnew = c(gcsnew,mean(gccurr[,2],na.rm=T))
  }
  gcinfo_grouped = data.frame(V1 = groups, V2 = gcsnew)
  write.table(gcinfo_grouped,gcpathnew,quote=F,sep="\t",row.names=F,col.names=F)
  #### group the coverage file:
  tssinfoname = basename(tssinfopath)
  nn = gsub(".txt","",tssinfoname)
  file2write = paste(outdir,"/",nn,".overlaping.bases.info.file.txt",sep="")
  overlapnewpath = gsub(".txt",paste(".merged.by.",gsize,".txt",sep=""),file2write)
  overlapinfo = fread(file2write,data.table=F)
  groups = as.matrix(grouped_info$geneName)
  SelectorNDRbasesall = c()
  Selector2Kbasesall = c()
  chrs = c()
  starts = c(); ends = c(); tss_ids = c(); genenames = c()
  for (g in groups)
  {
    geneshere = unlist(strsplit(as.character(g),"_&&_"))
    ocurr = overlapinfo[match(geneshere,overlapinfo$genename),]
    SelectorNDRbasesall = c(SelectorNDRbasesall,sum(ocurr$SelectorNDRbases))
    Selector2Kbasesall = c(Selector2Kbasesall,sum(ocurr$Selector2Kbases))
    chrs = c(chrs, paste(ocurr$chr,collapse="_"))
    starts = c(starts, paste(ocurr$start,collapse="_"))
    ends = c(ends, paste(ocurr$end,collapse="_"))
    genenames = c(genenames,paste(ocurr$genename,collapse="_&&_"))
    tss_ids = c(tss_ids, paste(ocurr$genename,collapse="_&&_"))
  }
  overlap_info = data.frame(ch = chrs, start = starts, end = ends, genename = genenames,
                            SelectorNDRbases = SelectorNDRbasesall, Selector2Kbases = Selector2Kbasesall,
                            tss_id = tss_ids)
  write.table(overlap_info,overlapnewpath,quote=F,sep="\t",row.names=F)
  return(list(gepref = grouped_info, gcpath = gcpathnew, agg_gc = gcinfo_grouped, overlapfile = overlapnewpath))
}