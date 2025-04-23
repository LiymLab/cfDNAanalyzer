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

priordatadir =paste0(path2dep,"/priordata/")
##################################################################################
rnaseqpath = paste0(priordatadir,"/reference.pbmc.RNA.TPM.txt")
rnaseqpath2 = paste0(priordatadir,"/reference.pbmc.RNA.TPM.ExomeSpace.txt")
gcpath = paste0(priordatadir,"/ensembl75.TSS_GC_aggregate.final.txt")
##################################################################################
pbmc_genome = fread(rnaseqpath,data.table = F)
pbmc_exome = fread(rnaseqpath2,data.table = F)
######################################################################
gc_file = read.table(gcpath,header=T,sep="\t")
agg_gc = gc_file
#####################################################################
generate_insilico <- function(mat_in, nfrag)
{
  mat_in_new = mat_in
  allsamples = apply(mat_in,1,function(x){z=x/sum(x);return(rmultinom(n = 1, prob = z,size =nfrag ))})
  mat_in_new = (t(allsamples))
  colnames(mat_in_new) = colnames(mat_in)
  return(mat_in_new)
}
#####################################################################
givemematrix <- function(df_in, gepref)
{
  Xin = t(t(df_in[, paste("bin",seq(50,390),sep="_")]))  # 
  myintersected_genes = intersect(df_in$genename,gepref$gene)
  gepref_gtex_intersect = gepref[match(myintersected_genes,gepref$gene),]
  #######################################################################
  gepref_gtex_intersect_ordered = gepref_gtex_intersect[order(-gepref_gtex_intersect$avg_exp),]
  sample_sorted = Xin[match(gepref_gtex_intersect_ordered$gene,df_in$genename),]
  deep_wgs_sorted = df_in[match(gepref_gtex_intersect_ordered$gene,df_in$genename),]
  mynewsummary = givemestats(sample_sorted,gepref_gtex_intersect_ordered, groupsize=1, logspace = FALSE, order = TRUE)
  mynewX = mynewsummary$grouped_matrix
  mysamplemode = seq(50,390)[which.max(apply(mynewX,2,sum))]
  mynewX_adjusted = as.data.frame(mynewX)
  geneinblood = gepref_gtex_intersect_ordered$gene[order(-gepref_gtex_intersect_ordered[,2])]
  return(list(mat = mynewX_adjusted, genes = geneinblood, exp = mynewsummary$avg_exp))
}
########################################################################
allbins = seq(50,398)
allbins_str = paste("bin_",allbins, sep="")
entbins = seq(100,300)
entbins_str = paste("bin_",entbins, sep="")
##############
bins_to150 = seq(65,150)
bins_to150_str = paste("bin_",bins_to150, sep="")
bins_150to200 = seq(151,200)
bins_150to200_str = paste("bin_",bins_150to200, sep="")
bins_201to = seq(201,398)
bins_201to_str = paste("bin_",bins_201to, sep="")
window = 5
print("Initialization is done")

