##----------------------------------------------------------------
## Copyright Stanford University 2021
## M.S.Esfahani, et al., Alizadeh and Diehn labs, Stanford University
##----------------------------------------------------------------
epic_calcocf <- function(allsigs,gsize=1)
{
  mysignal_features_all = c()
  sigtype_cnt = 0
  for (sigpathname in names(allsigs))
  {
    sigpaths = allsigs[[sigpathname]]
    savestrname = sigpathname
    sigtype_cnt = sigtype_cnt + 1
    sgpth = sigpaths[1]
    signal_wgsHere = fread(sgpth,data.table = F)
    if (gsize<=1)
    {
      signal_wgsHere$start = signal_wgsHere$TSS-1000
      signal_wgsHere$end = signal_wgsHere$TSS+1000
      final_dataSig = renameTSS(signal_wgsHere)
    }else{
      u = as.numeric(unlist(lapply(signal_wgsHere$TSS,function(x){z=as.character(x);u=unlist(strsplit(z,"_"))[1];return(u)})))
      signal_wgsHere$start = u -1000
      signal_wgsHere$end = u +1000
      final_dataSig = signal_wgsHere
      final_dataSig$tss_id = signal_wgsHere$genename
    }
    genesHere = final_dataSig$genename
    tssids = final_dataSig$tss_id
    mysignal_features = data.frame(gene = genesHere, TSS_ID = tssids, gc = NA,
                                   OCF_x = final_dataSig$OCF_x, 
                                   OCF_1 = final_dataSig$OCF_1,
                                   OCF_LDiff = final_dataSig$OCF_LDiff,
                                   OCF_LTot = final_dataSig$OCF_LTot,
                                   OCF_RDiff = final_dataSig$OCF_RDiff,
                                   OCF_RTot = final_dataSig$OCF_RTot)
    mysignal_features$gc = agg_gc$gc[match(mysignal_features$gene,agg_gc$gene)]
    colnames(mysignal_features)[c(4,5,6,7,8,9)] = 
      paste(sigpathname,colnames(mysignal_features)[c(4,5,6,7,8,9)],sep="_")
    if (sigtype_cnt==1)
    {
      mysignal_features_all = mysignal_features
    }else{
      mysignal_features_all = merge(mysignal_features_all,mysignal_features,by=c("TSS_ID","gc","gene"))
    }
  }
  return(mysignal_features_all)
}
