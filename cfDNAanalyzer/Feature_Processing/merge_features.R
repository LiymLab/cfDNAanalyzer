#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
input_dir = args[1]
output_dir = args[2]
Features = args[3]

if ("EM" %in% strsplit(Features, ",")[[1]]) {
  EM_motifs_frequency_path = paste0(input_dir, "/EM_motifs_frequency.csv")
  if (!file.exists(EM_motifs_frequency_path)) {
    warning(paste("File not found:", EM_motifs_frequency_path))
    return(NULL)
  }
  
  EM_motifs_mds_path = paste0(input_dir, "/EM_motifs_mds.csv")
  if (!file.exists(EM_motifs_mds_path)) {
    warning(paste("File not found:", EM_motifs_mds_path))
    return(NULL)
  }
  
  
  EM_motifs_frequency = read.csv(EM_motifs_frequency_path, header = T)
  EM_motifs_mds = read.csv(EM_motifs_mds, header = T)
  EM = merge(EM_motifs_frequency, EM_motifs_mds, by=c("sample","label"))
  write.csv(EM, file= paste0(output_dir, "/EM.csv"), quote = F, row.names = F)
}


if ("EMR" %in% strsplit(Features, ",")[[1]]) {
  EMR_aggregated_mds_path = paste0(input_dir, "/EMR_aggregated_mds.csv")
  if (!file.exists(EMR_aggregated_mds_path)) {
    warning(paste("File not found:", EMR_aggregated_mds_path))
    return(NULL)
  }
  
  EMR_aggregated_motif_frequency_path = paste0(input_dir, "/EMR_aggregated_motif_frequency.csv")
  if (!file.exists(EMR_aggregated_motif_frequency_path)) {
    warning(paste("File not found:", EMR_aggregated_motif_frequency_path))
    return(NULL)
  }
  
  EMR_region_mds_path = paste0(input_dir, "/EMR_region_mds.csv")
  if (!file.exists(EMR_region_mds_path)) {
    warning(paste("File not found:", EMR_region_mds_path))
    return(NULL)
  }
  
  EMR_region_motif_frequency_path = paste0(input_dir, "/EMR_region_motif_frequency.csv")
  if (!file.exists(EMR_region_motif_frequency_path)) {
    warning(paste("File not found:", EMR_region_motif_frequency_path))
    return(NULL)
  }

  
  EMR_aggregated_mds = read.csv(EMR_aggregated_mds_path, header = T)
  EMR_aggregated_motif_frequency = read.csv(EMR_aggregated_motif_frequency_path, header = T)
  EMR_region_mds = read.csv(EMR_region_mds_path, header = T)
  EMR_region_motif_frequency = read.csv(EMR_region_motif_frequency_path, header = T)
  
  EMR = merge(EMR_aggregated_mds, EMR_aggregated_motif_frequency, by=c("sample","label"))
  EMR = merge(EMR, EMR_region_mds, by=c("sample","label"))
  EMR = merge(EMR, EMR_region_motif_frequency, by=c("sample","label"))
  write.csv(EMR, file= paste0(output_dir, "/EMR.csv"), quote = F, row.names = F)
}


if ("NOF" %in% strsplit(Features, ",")[[1]]) {
  NOF_meanfuziness_path = paste0(input_dir, "/NOF_meanfuziness.csv")
  if (!file.exists(NOF_meanfuziness_path)) {
    warning(paste("File not found:", NOF_meanfuziness_path))
    return(NULL)
  }
  
  NOF_occupancy_path = paste0(input_dir, "/NOF_occupancy.csv")
  if (!file.exists(NOF_occupancy_path)) {
    warning(paste("File not found:", NOF_occupancy_path))
    return(NULL)
  }
  
  
  NOF_meanfuziness = read.csv(NOF_meanfuziness_path, header = T)
  colnames(NOF_meanfuziness)[-c(1, 2)] <- paste0(colnames(NOF_meanfuziness)[-c(1, 2)], "_meanfuziness")
  NOF_occupancy = read.csv(NOF_occupancy_path, header = T)
  colnames(NOF_occupancy)[-c(1, 2)] <- paste0(colnames(NOF_occupancy)[-c(1, 2)], "_occupancy")
  
  NOF = merge(NOF_meanfuziness, NOF_occupancy, by=c("sample","label"))
  write.csv(NOF, file= paste0(output_dir, "/NOF.csv"), quote = F, row.names = F)
}


if ("NP" %in% strsplit(Features, ",")[[1]]) {
  NP_amplitude_path = paste0(input_dir, "/NP_amplitude.csv")
  if (!file.exists(NP_amplitude_path)) {
    warning(paste("File not found:", NP_amplitude_path))
    return(NULL)
  }
  
  NP_central_coverage_path = paste0(input_dir, "/NP_central_coverage.csv")
  if (!file.exists(NP_central_coverage_path)) {
    warning(paste("File not found:", NP_central_coverage_path))
    return(NULL)
  }
  
  NP_mean_coverage_path = paste0(input_dir, "/NP_mean_coverage.csv")
  if (!file.exists(NP_mean_coverage_path)) {
    warning(paste("File not found:", NP_mean_coverage_path))
    return(NULL)
  }
  
  NP_amplitude = read.csv(NP_amplitude_path, header = T)
  colnames(NP_amplitude)[-c(1, 2)] <- paste0(colnames(NP_amplitude)[-c(1, 2)], "_amplitude")
  NP_central_coverage = read.csv(NP_central_coverage_path, header = T)
  colnames(NP_central_coverage)[-c(1, 2)] <- paste0(colnames(NP_central_coverage)[-c(1, 2)], "_central_coverage")
  NP_mean_coverage = read.csv(NP_mean_coverage_path, header = T)
  colnames(NP_mean_coverage)[-c(1, 2)] <- paste0(colnames(NP_mean_coverage)[-c(1, 2)], "_mean_coverage")
  
  NP = merge(NP_amplitude, NP_central_coverage, by=c("sample","label"))
  NP = merge(NP, NP_mean_coverage, by=c("sample","label"))
  write.csv(NP, file= paste0(output_dir, "/NP.csv"), quote = F, row.names = F)
}


if ("WPS" %in% strsplit(Features, ",")[[1]]) {
  WPS_long_path = paste0(input_dir, "/WPS_long.csv")
  if (!file.exists(WPS_long_path)) {
    warning(paste("File not found:", WPS_long_path))
    return(NULL)
  }
  
  WPS_short_path = paste0(input_dir, "/WPS_short.csv")
  if (!file.exists(WPS_short_path)) {
    warning(paste("File not found:", WPS_short_path))
    return(NULL)
  }
  
  
  WPS_long = read.csv(WPS_long_path, header = T)
  colnames(WPS_long)[-c(1, 2)] <- paste0(colnames(WPS_long)[-c(1, 2)], "_long")
  WPS_short = read.csv(WPS_short_path, header = T)
  colnames(WPS_short)[-c(1, 2)] <- paste0(colnames(WPS_short)[-c(1, 2)], "_short")
  
  WPS = merge(WPS_long, WPS_short, by=c("sample","label"))
  write.csv(WPS, file= paste0(output_dir, "/WPS.csv"), quote = F, row.names = F)
}
