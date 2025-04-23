##----------------------------------------------------------------
## Copyright Stanford University 2021
## M.S.Esfahani, et al., Alizadeh and Diehn labs, Stanford University
##----------------------------------------------------------------

epic_calcndr <- function(allsigs,gsize = 1)
{
  allbins = seq(50,398)
  allbins_str = paste("bin_",allbins, sep="")
  allbins4ndr = seq(min(entbins),max(entbins))
  allbins4ndr_str = paste("bin_",allbins4ndr, sep="")
  bins_to150 = seq(65,150)
  bins_to150_str = paste("bin_",bins_to150, sep="")
  bins_150to200 = seq(151,200)
  bins_150to200_str = paste("bin_",bins_150to200, sep="")
  bins_201to = seq(201,398)
  bins_201to_str = paste("bin_",bins_201to, sep="")
  mysignal_features_all = c()
  sigtype_cnt = 0
  aggby = 5
  for (sigpathname in names(allsigs))
  {
    sigpaths = allsigs[[sigpathname]]
    savestrname = sigpathname
    sigtype_cnt = sigtype_cnt + 1
    sgcn = 0
    for (sgpth in sigpaths)
    {
      signal_wgsHere = fread(sgpth,data.table = F)
      nnames_here = paste(signal_wgsHere$genename,signal_wgsHere$chr,signal_wgsHere$start,sep="_")
      sgcn=sgcn+1
      if (sgcn==1)
      {
        signal_wgsFinal = signal_wgsHere
      }else{
        nnames_final = paste(signal_wgsFinal$genename,signal_wgsFinal$chr,signal_wgsFinal$start,sep="_")
        nnames_intersect = intersect(nnames_final,nnames_here)
        signal_wgsFinal[allbins_str] =  (signal_wgsFinal[allbins_str])+ (signal_wgsHere[allbins_str])
      }
    }
    ##################################### SUBTRACTING THE FLANKING REGIONS
    signal_wgsFinal[allbins_str] = signal_wgsFinal[allbins_str]
    #####################################
    signal_wgs = signal_wgsFinal
    signal_wgs = signal_wgs
    nnames_sig = paste(signal_wgs$genename,signal_wgs$chr,signal_wgs$start,sep="_")
    signal_wgs = signal_wgs
    newsignal = signal_wgs
    samplenames = c()
    count_here = 0
    final_dataSig = signal_wgs
    if (gsize<=1)
    {
      final_dataSig = renameTSS(final_dataSig)
    }else{
      final_dataSig$tss_id = final_dataSig$genename
      
    }
    all_curr_matrix = final_dataSig[allbins_str]
    genesHere = final_dataSig$genename
    tssids = final_dataSig$tss_id
    ### all bins
    allbins_sum = apply(all_curr_matrix[allbins4ndr_str],2,sum,na.rm=T)
    fulldepth = sum(allbins_sum)
    totDepthSig = apply(all_curr_matrix[allbins4ndr_str],1,sum,na.rm=T)
    totDepthSigNorm_allbins = totDepthSig
    ### bins_to150_str
    allbins_sum = apply(all_curr_matrix[bins_to150_str],2,sum,na.rm=T)
    fulldepth = sum(allbins_sum)
    alldistsHERE = t(apply(all_curr_matrix[bins_to150_str],1,convert2dist,ww = window))
    totDepthSig = apply(all_curr_matrix[bins_to150_str],1,sum,na.rm=T)
    totDepthSigNorm_short = totDepthSig
    allbins_sum = apply(all_curr_matrix[bins_150to200_str],2,sum,na.rm=T)
    fulldepth = sum(allbins_sum)
    alldistsHERE = t(apply(all_curr_matrix[bins_150to200_str],1,convert2dist,ww = window))
    totDepthSig = apply(all_curr_matrix[bins_150to200_str],1,sum,na.rm=T)
    totDepthSigNorm_middle = totDepthSig
    ### bins_201to_str
    allbins_sum = apply(all_curr_matrix[bins_201to_str],2,sum,na.rm=T)
    fulldepth = sum(allbins_sum)
    alldistsHERE = t(apply(all_curr_matrix[bins_201to_str],1,convert2dist,ww = window))
    totDepthSig = apply(all_curr_matrix[bins_201to_str],1,sum,na.rm=T)
    totDepthSigNorm_long = totDepthSig
    mediandepth = median(totDepthSig)
    totreads = sum(apply(signal_wgs[allbins_str],2,sum))
    mysignal_features = data.frame(gene = genesHere, TSS_ID = tssids, 
                                   gc = NA,                                   
                                   NDR_all= totDepthSigNorm_allbins, 
                                   NDR_short = totDepthSigNorm_short,
                                   NDR_middle = totDepthSigNorm_middle,
                                   NDR_long = totDepthSigNorm_long)
    
    mysignal_features$gc = agg_gc$gc[match(mysignal_features$gene,agg_gc$gene)]
    mysignal_features$Total_reads = totreads
    colnames(mysignal_features)[c(4,5,6,7,8)] = paste(sigpathname,colnames(mysignal_features)[c(4,5,6,7,8)],sep="_")
    if (sigtype_cnt==1)
    {
      mysignal_features_all = mysignal_features
    }else{
      mysignal_features_all = merge(mysignal_features_all,mysignal_features,by=c("TSS_ID","gc","gene"))
    }
    
  }
  return(mysignal_features_all)
}

