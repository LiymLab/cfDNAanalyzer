##----------------------------------------------------------------
## Copyright Stanford University 2021
## M.S.Esfahani, et al., Alizadeh and Diehn labs, Stanford University
##----------------------------------------------------------------



samplewise_adjust_final <- function(dfin,neggenes){
  negent = median(dfin$Shannon[dfin$gene%in%neggenes])
  xtop = quantile(dfin$Shannon,0.95)
  xbottom = quantile(dfin$Shannon,0.05)
  scalingfactor = xtop/negent
  currmean = median(dfin$NDR_relative_2k[dfin$gene%in%neggenes])
  dfin$PFE = pmin(1,dfin$PFE/pgamma(scalingfactor-1,shape= 0.5,rate = 1))
  dfin$NDR_relative_2k = (dfin$NDR_relative_2k/currmean)
  return(dfin)
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
renameTSS <- function(df)
{
  uniqgenes = unique(df$genename)
  numtss = table(df$genename)
  singletons = which(numtss==1)
  nonsingle = which(numtss>1)
  df$tss_id = df$genename
  df$tss_id[singletons] = paste(df$genename[singletons],1,sep="_")
  
  if (length(nonsingle)>0)
  {
    
    nonsingleidx = which(df$genename%in%names(numtss)[nonsingle])
    allgenenames = df$genename
    allstart = df$start
    uqnonsingle = unique(df$genename[nonsingleidx])
    for (g in uqnonsingle)
    {
      idxcur = which(allgenenames==g)
      df$tss_id[idxcur] = paste(df$genename[idxcur],seq(length(idxcur))[order(allstart[idxcur])],sep="_")
    }
  }
  return(df)
}

creatergnfile <- function(tsspath,uu,dd,outdir)
{
  # 20241105 zjp modified
  options(scipen = 200)
  # end
  library(data.table)
  alltss = fread(tsspath, data.table = F)
  colnames(alltss) = c("CHR","TSS","Gene-Symbol","Category","Strand","TSS_ID")
  nnn = nrow(alltss)
  alltss$Start = alltss$TSS - (uu*(alltss$Strand==+1)+dd*(alltss$Strand==-1))
  alltss$End = alltss$TSS + (dd*(alltss$Strand==+1)+uu*(alltss$Strand==-1))
  filename = basename(tsspath)
  tsspath_new = gsub(".txt",paste("_upstream_",uu,"_downstream_",dd,".txt",sep=""),paste(outdir,"/",filename,sep=""))
  #
  write.table(alltss[c("CHR","Start","End","Gene-Symbol","Category","TSS")], tsspath_new ,col.names = F,
              quote=F,sep="\t",row.names=F)
  return(tsspath_new)
}


rmvgenes <- function(dfin, path2file = 'DLBCL_grant/tss_rgns_allgenes_exclude.uniq.txt')
{
  ## genes to remove
  gene2rmv_df = fread(path2file,data.table=F)
  gene2rmv = unique(gene2rmv_df$V3)
  dfin = dfin[!dfin$genename%in%gene2rmv,]
  return(dfin)
}
gc.correct <- function(coverage, bias, every = 10) {
  coverage_sorted = coverage[order(-bias)]
  bias_sorted = bias[order(-bias)]
  coverage_sorted_smooth = rollapply(coverage_sorted,every,by=every,mean)
  bias_sorted_smooth = rollapply(bias_sorted,every,by=every,mean)
  mycor1 = cor.test(coverage_sorted_smooth,bias_sorted_smooth, method = "spearman")
  mycor2 = cor.test(coverage_sorted_smooth,bias_sorted_smooth, method = "pearson")
  if (!is.na(mycor1) & !is.na(mycor2))
  {
    if (max(mycor1$estimate,mycor2$estimate)>0.2)
      # 
    {
      i <- seq(min(bias_sorted_smooth, na.rm=TRUE), max(bias_sorted_smooth, na.rm=TRUE), by = 0.005)
      coverage.trend <- loess(coverage_sorted_smooth ~ bias_sorted_smooth,span=0.5)
      coverage.model <- loess(predict(coverage.trend, i) ~ i,span=0.50)
      coverage.pred <- predict(coverage.model, bias)
      coverage.corrected <- coverage - coverage.pred + median(coverage,na.rm=T)
      coverage.corrected[is.na(coverage.corrected)] = median(coverage,na.rm=T)
    }else{
      coverage.corrected = coverage
    }
  }else{
    coverage.corrected = coverage
    
  }
  return(coverage.corrected)
}


keep1TSS <- function(dfin, entbins = seq(65,320))
{
  strcols = paste("bin_",entbins,sep="")
  depths = apply(dfin[strcols],1,sum,na.rm=T)
  overallsum =  apply(dfin[strcols],2,sum,na.rm=T)
  myfulldist = convert2dist(overallsum,ww = 3)
  quickEnt = apply(dfin[strcols],1,function(x){t = calculate_ent_bayesian_avg(x,distfull =myfulldist,
                                                                              alpha0 =median(depths,na.rm=T));return(t)})
  quickEnt[is.infinite(quickEnt)] = 1
  quickEnt[is.na(quickEnt)] = 1
  dfin$geneEnt = paste(dfin$genename,quickEnt,sep="_")
  agg = aggregate(quickEnt, by = list(genename= dfin$genename), max)
  agginfo = paste(agg$genename,agg$x,sep="_")
  dfin_new = (dfin[match(agginfo,dfin$geneEnt),])
  return(dfin_new)
}
fim <- function(x){
  xx = sum((diff(sqrt(x)))^2)
  return(xx)
}
renyi <- function(x,alpha=2){
  xx = 1/(1-alpha)*log2(sum(x^alpha))
  return(xx)
}
gini_index <- function(x,range = seq(100,300)){
  xx = Gini(range,n=x)
  return(xx)
}
calculate_sd <- function(x,range){mu=sum(x*range);mom2 = sum(x*range^2);return(sqrt(mom2-mu^2))}
calculate_mu <- function(x,range){mu=sum(x*range);return(mu)}
calculate_ent <- function(x,range){x=x/sum(x);ent=-sum(x*log2(x));return(ent)}

calculate_ent_bayesian <- function(counts,distfull, alpha0,n=250, baseents = NULL){
  countsnew = counts
  if(sum(is.na(countsnew))>0)
  {
    print("calculate_ent_bayesian(): Encountered NA counts, returning PFE as NaN.")
    return(NaN)
  }
  distfullnew = distfull
  if(sum(is.na(distfullnew))>0)
  {
    print("calculate_ent_bayesian(): Encountered NA distfull, returning PFE as NaN.")
    return(NaN)
  }
  distfullnew = distfullnew/sum(distfullnew)
  updatedalpha = 1+(distfullnew*alpha0+countsnew)
  if(sum(is.na(updatedalpha))>0)
  {
    print("calculate_ent_bayesian(): Encountered NA alpha, returning PFE as NaN.")
    return(NaN)
  }
  if (sum(updatedalpha)<10)
  {
    updatedalpha = updatedalpha/sum(updatedalpha)*10
  }
  N = sum(counts)
  mydirichlet <- rdirichlet(n=n, alpha=as.numeric(updatedalpha)) + 1e-5
  mynewent = digamma(sum(updatedalpha)+1) - sum(updatedalpha/sum(updatedalpha)*digamma(updatedalpha+1))
  newnet = mynewent/log(2)
  base_ent = calculate_ent(distfullnew)
  gammanew <- function(arrayin){o=pgamma(arrayin,shape=.5,rate=1);return(o)}
  if (is.null(baseents))
  {
    allents = (apply(mydirichlet,1,function(x){x=x+1e-5;
    ent = calculate_ent(x/sum(x));return((ent/base_ent))}))- 1
    bayesian_sig = mean(gammanew(allents))
  }else{
    bayesian_sig = 0
    dirichletEnts = (apply(mydirichlet,1,function(x){x=x+1e-5;ent = calculate_ent(x/sum(x));return((ent))}))
    for (base_ent in baseents)
    {
      allents = dirichletEnts/base_ent - 1
      bayesian_sig = bayesian_sig+mean(gammanew(allents))
    }
    bayesian_sig = bayesian_sig/length(baseents)
  }
  return(bayesian_sig)
}
# calculate_ent_bayesian <- function(counts,distfull, alpha0,n=250, baseents = NULL){
#   countsnew = counts
#   distfullnew = distfull
#   distfullnew = distfullnew/sum(distfullnew)
#   updatedalpha = 1+(distfullnew*alpha0+countsnew)
#   if (sum(updatedalpha)<10)
#   {
#     updatedalpha = updatedalpha/sum(updatedalpha)*10
#   }
#   N = sum(counts)
#   mydirichlet <- rdirichlet(n=n, alpha=as.numeric(updatedalpha)) + 1e-5
#   mynewent = digamma(sum(updatedalpha)+1) - sum(updatedalpha/sum(updatedalpha)*digamma(updatedalpha+1))
#   newnet = mynewent/log(2)
#   base_ent = calculate_ent(distfullnew)
#   gammanew <- function(arrayin){o=pgamma(arrayin,shape=.5,rate=1);return(o)}
#   if (is.null(baseents))
#   {
#     allents = (apply(mydirichlet,1,function(x){x=x+1e-5;
#     ent = calculate_ent(x/sum(x));return((ent/base_ent))}))- 1
#     bayesian_sig = mean(gammanew(allents))
#   }else{
#     bayesian_sig = 0
#     dirichletEnts = (apply(mydirichlet,1,function(x){x=x+1e-5;ent = calculate_ent(x/sum(x));return((ent))}))
#     for (base_ent in baseents)
#     {
#       allents = dirichletEnts/base_ent - 1
#       bayesian_sig = bayesian_sig+mean(gammanew(allents))
#     }
#     bayesian_sig = bayesian_sig/length(baseents)
#   }
#   return(bayesian_sig)
# }


myrollapply <- function(u,ww){rollapply(u,ww,mean,fill=0)}

convert2dist <- function(y,ww = 5)
{
  y=c(y,rep(y[length(y)],ww-1));
  myz=(y+0.0001)/sum(y+0.0001)
  myz = myrollapply(myz,ww=ww)
  myz[myz<=1e-7] = 1e-7
  return(myz)
}
givemestats <- function(initialX_reordered,gepref, groupsize=10, logspace = FALSE,
                        psudoadd = 0.1, order = TRUE)
{
  
  X_grouped = c()
  avg_exp = c()
  newgenenames = c()
  if (order)
  {
    initialX_reordered = initialX_reordered[order(-gepref[,2]),]
    gepref_ordered = gepref[order(-gepref[,2]),]
    avg_exp = gepref_ordered[,2]
  }else{
    initialX_reordered = initialX_reordered
    gepref_ordered = gepref
    avg_exp = gepref_ordered[,2]
    newgenenames = gepref_ordered[,1]
  }
  rawgenes = gepref[,1]
  numg = groupsize
  if (numg>1)
  {
    avg_exp = c()
    newgenenames = c()
    for (kk in seq(1,floor(length(rawgenes)/numg)))
    {
      myspace = seq((kk-1)*numg+1,kk*numg)
      X_sub = colSums(initialX_reordered[myspace,])
      X_grouped = rbind(X_grouped, X_sub)
      avg_exp = c(avg_exp,mean(gepref_ordered[myspace,2]))
      newgenenames = c(newgenenames, paste(gepref_ordered[myspace,1],collapse="_"))
    }
  }else
  {
    X_grouped = initialX_reordered
    newgenenames = gepref_ordered[,1]
  }
  X_grouped_normalized = t(apply(X_grouped,1,function(x){return((x+psudoadd)/sum(x+psudoadd))}))
  if (logspace)
  {
    mystats = list(means = colMeans(log2(X_grouped_normalized)), sds = apply(log2(X_grouped_normalized),2,sd),
                   grouped_matrix = X_grouped, avg_exp = avg_exp)
  }else{
    mystats = list(means = colMeans(X_grouped_normalized), sds = apply(X_grouped_normalized,2,sd), 
                   grouped_matrix = X_grouped, avg_exp = avg_exp, genenames = newgenenames)
    
  }
  return(mystats)
}

calculate_ent_bayesian_avg <- function(counts,distfull, alpha0){
  countsnew = counts
  if(sum(is.na(countsnew))>0)
  {
    print("calculate_ent_bayesian_avg(): Encountered NA counts, returning Shannon as NaN.")
    return(NaN)
  }
  distfullnew = distfull
  if(sum(is.na(distfullnew))>0)
  {
    print("calculate_ent_bayesian_avg(): Encountered NA distfull, returning Shannon as NaN.")
    return(NaN)
  }
  n = 250
  distfullnew = distfullnew/sum(distfullnew)
  updatedalpha = (distfullnew*alpha0+countsnew)
  if(sum(is.na(updatedalpha))>0)
  {
    print("calculate_ent_bayesian_avg(): Encountered NA alhpa, returning Shannon as NaN.")
    return(NaN)
  }
  if (sum(updatedalpha)<10)
  {
    updatedalpha = updatedalpha/sum(updatedalpha)*10
  }
  N = sum(counts)
  mynewent = digamma(sum(updatedalpha)+1) - sum(updatedalpha/sum(updatedalpha)*digamma(updatedalpha+1))
  newnet = mynewent/log(2)
  return(newnet)
}
# calculate_ent_bayesian_avg <- function(counts,distfull, alpha0){
#   countsnew = counts
#   distfullnew = distfull
#   n = 250
#   distfullnew = distfullnew/sum(distfullnew)
#   updatedalpha = (distfullnew*alpha0+countsnew)
#   if (sum(updatedalpha)<10)
#   {
#     updatedalpha = updatedalpha/sum(updatedalpha)*10
#   }
#   N = sum(counts)
#   mynewent = digamma(sum(updatedalpha)+1) - sum(updatedalpha/sum(updatedalpha)*digamma(updatedalpha+1))
#   newnet = mynewent/log(2)
#   return(newnet)
# }
