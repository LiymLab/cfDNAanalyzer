# Usage:
# 1. conda activate cfDNATools_new
# 2. Rscript install_R_packages.R

# IMPORTANT: You may want to choose a CRAN mirror that is close to your geographical location for faster downloads.
# This script uses the Tsinghua University mirror in China, which is fast for users in China.
# You can find a list of CRAN mirrors at: https://cran.r-project.org/mirrors.html
# and replace the URL below with your preferred mirror.

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))


if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# CRAN packages list
cran_packages_versions <- list(
  DescTools = "0.99.40",
  zoo = "1.8.12",
  plyr = "1.8.9",
  reshape2 = "1.4.4",
  data.table = "1.15.2",
  MASS = "7.3-60.0.1",
  e1071 = "1.7-14",
  gtools = "3.9.5",
  matrixStats = "1.2.0",
  optparse = "1.7.4",
  httr = "1.4.7",
  tidyverse = "2.0.0",
  RCurl = "1.98-1.14"
)

for (pkg in names(cran_packages_versions)) {
  remotes::install_version(pkg, version = cran_packages_versions[[pkg]])
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(version = "3.18")

# Bioconductor packages list
bioconductor_packages <- c("HMMcopy", "GenomeInfoDb", "GenomicRanges", 
                           "Rsamtools", "GenomicAlignments", "biovizBase", 
                           "BSgenome.Hsapiens.UCSC.hg19","BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install(bioconductor_packages)

installed <- rownames(installed.packages())
missing_packages <- setdiff(c(names(cran_packages_versions), bioconductor_packages), installed)

if (length(missing_packages) > 0) {
  cat("The following packages were not successfully installed:", paste(missing_packages, collapse = ", "), "\n")
} else {
  cat("All packages were successfully installed.\n")
}
