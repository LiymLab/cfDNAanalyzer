#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
input_dir = args[1]
output_dir = args[2]

label_file = read.table(paste0(input_dir,"/label.txt"), header = T, sep = ",")

all_dirs <- list.dirs(path = input_dir , full.names = FALSE, recursive = FALSE)
sub_dirs <- all_dirs[all_dirs != "."]
folder_names <- sub_dirs

for (i in c(1:377)) {
    all_data = data.frame()
    
    for (j in folder_names) {
        cov_file <- file.path(input_dir, sample_dir, "NP", "sample_name_1", 
                              "sample_name_1.GC_corrected.coverage.txt")

        if (!file.exists(cov_file)) {
          warning(paste("File not found:", cov_file))
          return(NULL)
        }
        
        data <-  read.table(cov_file, header = TRUE)
        
        if (nrow(data) < 377) {
          warning(paste("Not 377 rows in", cov_file))
          return(NULL)
        }
        
        data = data[i,]
        label <- label_file[label_file$sample %in% j, "label"]
        data$sample = j
        data$label = label
        
        data <- data[, c("sample", "label", setdiff(names(data), c("sample", "label")))]
        all_data = rbind(all_data, data)
    }
    
    write.table(all_data, file = paste0(output_dir, "/site_list_", i, ".txt"), quote = F, sep = ",", row.names = F)
}

