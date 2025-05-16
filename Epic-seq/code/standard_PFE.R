#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
input_file = args[1]
output_file = args[2]

data = read.csv(input_file)
sample_col = data$sample
label_col = data$label
data$sample = NULL
data$label = NULL
df_row_scaled <- t(apply(data, 1, scale))
df_row_scaled = as.data.frame(df_row_scaled)
colnames(df_row_scaled) <- colnames(data)
df_row_scaled$sample = sample_col
df_row_scaled$label = label_col
df_row_scaled <- df_row_scaled[, c("sample", "label", setdiff(names(df_row_scaled), c("sample", "label")))]
write.csv(df_row_scaled, file = output_file, quote = F, row.names = F)
