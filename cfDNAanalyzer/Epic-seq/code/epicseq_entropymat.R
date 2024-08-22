##----------------------------------------------------------------
## Copyright Stanford University 2021
## M.S.Esfahani, et al., Alizadeh and Diehn labs, Stanford University
##----------------------------------------------------------------

epic_calcpfe <- function(allsigs,allbgs,gsize=1,mode="epicseq")
{
  mysignal_features_all = c()
  sigtype_cnt = 0
  for (sigpathname in names(allsigs))
  {
    sigpaths = allsigs[[sigpathname]]
    savestrname = sigpathname
    sigtype_cnt = sigtype_cnt + 1
    #####################################################################
    bgpaths = allbgs[[sigpathname]]
    #####################################################################
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
    #### working on background
    bgcn = 0
    for (bgpath in bgpaths)
    {
      bg_wgsHere = fread(bgpath,data.table = F)
      nnames_here = paste(bg_wgsHere$genename,bg_wgsHere$chr,bg_wgsHere$start,sep="_")
      bgcn=bgcn+1
      if (bgcn==1)
      {
        bg_wgsFinal = bg_wgsHere
      }else{
        nnames_final = paste(bg_wgsFinal$genename,bg_wgsFinal$chr,bg_wgsFinal$start,sep="_")
        nnames_intersect = intersect(nnames_final,nnames_here)
        bg_wgsFinal[allbins_str] =  (bg_wgsFinal[allbins_str])+ (bg_wgsHere[allbins_str])
      }
    }
    ##################################### SUBTRACTING THE FLANKING REGIONS
    signal_wgsFinal[allbins_str] = signal_wgsFinal[allbins_str]
    signal_wgsFinal[allbins_str][signal_wgsFinal[allbins_str]<0] = 0
    #####################################
    signal_wgs = signal_wgsFinal
    newsignal = signal_wgs
    samplenames = c()
    count_here = 0
    final_dataSig = signal_wgs
    if (gsize<=1)
    {
      final_dataSig = renameTSS(final_dataSig)
      bg_wgsFinal = renameTSS(bg_wgsFinal)
      
    }else{
      final_dataSig$tss_id = final_dataSig$genename
      bg_wgsFinal$tss_id = bg_wgsFinal$genename
      
    }
    bg_wgsFinal_negCtr = bg_wgsFinal[grep("negativeControl",bg_wgsFinal$reason),]
    all_curr_matrix = final_dataSig[allbins_str]
    genesHere = final_dataSig$genename
    tssids = final_dataSig$tss_id
    allbins_sum = apply(all_curr_matrix[entbins_str],2,sum,na.rm=T)
    myfulldist = convert2dist(allbins_sum,ww = window)
    totDepthSig = apply(all_curr_matrix[entbins_str],1,sum,na.rm=T)
    totDepthSigNorm = totDepthSig/median(totDepthSig)
    totDepthSigNormLog = log2(0.001+totDepthSigNorm)
    totreads = sum(apply(signal_wgs[entbins_str],2,sum))
    ###### background distributions:
    if (mode=="epicseq")
    {
      nneg = min(20,nrow(bg_wgsFinal_negCtr))
    }else{
      nneg = nrow(bg_wgsFinal_negCtr)
    }
    set.seed(2020)
    allbins_sum_bg = colSums(bg_wgsFinal_negCtr[entbins_str],na.rm=T)
    alldists_bg = convert2dist(allbins_sum_bg, ww = window)
    allbins_sum_bg1 = colSums(bg_wgsFinal_negCtr[sample(seq(nneg),nneg,replace =T),entbins_str],na.rm=T)
    allbins_sum_bg2 = colSums(bg_wgsFinal_negCtr[sample(seq(nneg),nneg,replace =T),entbins_str],na.rm=T)
    allbins_sum_bg3 = colSums(bg_wgsFinal_negCtr[sample(seq(nneg),nneg,replace =T),entbins_str],na.rm=T)
    allbins_sum_bg4 = colSums(bg_wgsFinal_negCtr[sample(seq(nneg),nneg,replace =T),entbins_str],na.rm=T)
    
    baseents0 = calculate_ent(alldists_bg)
    baseents1 = calculate_ent(convert2dist(allbins_sum_bg1, ww = window))
    baseents2 = calculate_ent(convert2dist(allbins_sum_bg2, ww = window))
    baseents3 = calculate_ent(convert2dist(allbins_sum_bg3, ww = window))
    baseents4 = calculate_ent(convert2dist(allbins_sum_bg4, ww = window))
    baseents = c(baseents0,baseents1,baseents2,baseents3,baseents4)
    print("calculating PFE values...")
    print(dim(all_curr_matrix))
    bayeEntWGS <- function(x){calculate_ent_bayesian(x,distfull =alldists_bg[entbins_str],
                                                     alpha0 = 20, n = 200,
                                                     baseents= baseents)}
    
    bayeEntEPIC <- function(x){calculate_ent_bayesian(x,distfull =alldists_bg[entbins_str],
                                                      alpha0 =20, n = 300,
                                                      baseents= baseents)}
    ExpectedEntEPIC <- function(x){calculate_ent_bayesian_avg(x,distfull = myfulldist,
                                                              alpha0 =20)}
    if (mode =="wgs")
    {
      allentshere_adjustedSig = apply(all_curr_matrix[entbins_str],1,bayeEntWGS)
    }else{
      allentshere_adjustedSig = apply(all_curr_matrix[entbins_str],1,bayeEntEPIC)
    }
    shannonex = apply(all_curr_matrix[entbins_str],1,ExpectedEntEPIC)
    #allGini = apply(all_curr_matrix[entbins_str],1,function(x){y=gini_index(as.numeric(x),range = entbins);return(y)})
    # Chaorong: 20221103
    allGini = apply(all_curr_matrix[entbins_str],1,function(x){if(sum(is.na(x))>0){print("gini_index(): Encountered NA counts, returning GiniIndex as NaN.");return(NaN)};y=gini_index(as.numeric(x),range = entbins);return(y)})
    mysignal_features = data.frame(gene = genesHere, TSS_ID = tssids, gc = NA, 
                                   Depth= totDepthSigNormLog, 
                                   PFE = allentshere_adjustedSig,
                                   Shannon = shannonex,
                                   Gini = allGini)
    mysignal_features$gc = agg_gc$gc[match(mysignal_features$gene,agg_gc$gene)]
    mysignal_features$PFE = pmax(1e-5,mysignal_features$PFE)
    mysignal_features$PFE = pmin(1-1e-5,mysignal_features$PFE)
    colnames(mysignal_features)[c(4,5,6,7)] = paste(sigpathname,colnames(mysignal_features)[c(4,5,6,7)],sep="_")
    if (sigtype_cnt==1)
    {
      mysignal_features_all = mysignal_features
    }else{
      mysignal_features_all = merge(mysignal_features_all,mysignal_features,by=c("TSS_ID","gc","gene"))
    }
  }
  # 调试5
  print("Final signal features:")
  print(head(mysignal_features_all))
  #
  return(mysignal_features_all)
}