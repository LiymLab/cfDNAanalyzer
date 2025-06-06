##----------------------------------------------------------------
## Copyright Stanford University 2021
## M.S.Esfahani, et al., Alizadeh and Diehn labs, Stanford University
##----------------------------------------------------------------

group_by_gene <- function(fullfiles,str,type,gsize,outdir)
{
  allbins_str = paste("bin_",seq(50,398),sep="")
  gcpathfilename = basename(gcpath)
  gcpathnew = paste0(outdir,"/",gsub(".bed",paste(".grouped.by.",gsize,".bed",sep=""),gcpathfilename))  
  gcinfo = read.table(gcpathnew, header=F,sep="\t")
  groups = as.matrix(gcinfo[,1])
  ############# We now group all the histogram data---
  fnews = c()
  fcnt = 0
  for (f in fullfiles)
  {
    fcnt=fcnt+1
    xin = fread(f,data.table=F)
    fnew = gsub(".txt",paste(".grouped.by.",gsize,".txt",sep=""),f)
    allgeneshere = xin$genename
    dbins = c(); chrs = c(); starts = c(); ends = c();
    reasons = c(); strands = c(); samples = c(); 
    allvects = c()
    if (type!="ocf")
    {
      sample = xin$sampleName[1]
      if (fcnt==1)
      {
        for (g in groups)
        {
          geneshere = unlist(strsplit(as.character(g),"_&&_"))
          idxxx = which(allgeneshere%in%geneshere)
          xcurr = xin[idxxx,]
          xcurr = xcurr[match(geneshere,xcurr$genename),]
          vecs = colSums(xcurr[allbins_str],na.rm=T)
          chrs = c(chrs, paste(xcurr$chr,collapse="_"))
          starts = c(starts, paste(xcurr$start,collapse="_"))
          ends= c(ends,paste(xcurr$end,collapse="_"))
          reasons = c(reasons, paste(xcurr$reason,collapse="_"))
          strands = c(strands,paste(xcurr$strand,collapse="_"))
          allvects = rbind(allvects,vecs)
          
        }
      }else{
        for (g in groups)
        {
          geneshere = unlist(strsplit(as.character(g),"_&&_"))
          idxxx = which(allgeneshere%in%geneshere)
          xcurr = xin[idxxx,]
          xcurr = xcurr[match(geneshere,xcurr$genename),]
          vecs = colSums(xcurr[allbins_str],na.rm=T)
          allvects = rbind(allvects,vecs)
          
        }
      }
      
      allvects = as.data.frame(allvects)
      colnames(allvects) = allbins_str
      if (fcnt==1)
      {
        aug0 = data.frame(chr = chrs,start=starts, end = ends)
        aug1 = data.frame(genename = groups,reason=reasons, strand = strands)
      }
      aug1new = aug1
      aug1new$sampleName = rep(sample,nrow(aug1new))
      allnewdata = cbind(aug0,allvects,aug1new)
    }else{
      ocfcols = c("OCF_1","OCF_x", "OCF_LDiff","OCF_LTot","OCF_RDiff", "OCF_RTot")
      for (g in groups)
      {
        geneshere = unlist(strsplit(as.character(g),"_&&_"))
        xcurr = xin[match(geneshere,xin[,3]),]
        xcurrsum = xin[which(xin[,3]%in%geneshere),]
        chrs = c(chrs, paste(xcurr$chr,collapse="_"))
        starts = c(starts, paste(xcurr$TSS,collapse="_"))
        reasons = c(reasons, paste(xcurr[,12],collapse="_"))
        strands = c(strands,paste(xcurr[,4],collapse="_"))
        samples = c(samples, xcurr$sampleName[1])
        vecs = colSums(xcurrsum[ocfcols],na.rm=T)
        allvects = rbind(allvects,vecs)
      }
      allvects = as.data.frame(allvects)
      colnames(allvects) = ocfcols
      allnewdata = cbind(data.frame(chr = chrs,TSS=starts,genename = groups,strand=strands),allvects,
                         data.frame(reason=reasons, sampleName = samples ))
    }
    write.table(allnewdata,fnew,quote=F,sep="\t",row.names=F)
    fnews = c(fnews,fnew)
  }
  return(fnews)
}
