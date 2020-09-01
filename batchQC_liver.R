# run batchQC tests on liver data

library(BatchQC) # 4.0
library(devtools)
library(Biobase)
library(sva)
library(bladderbatch) # 4.0
library(snpStats) # 4.0
library(tximport)
library(DESeq2)
library(tidyverse)
library(ggfortify)
library(ggplot2)

## no negative values were changes to 0 before batchQC


## import counts tables
combat_counts <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Corrected/ComBat/Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
combat_seq_counts <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Corrected/ComBat_seq/Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
deseq2_counts <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Processed/All_Liver/DESeq2/DSeq2_NormCounts/Normalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
mbatch_eb_counts <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Corrected/MBatch_EB/Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
mbatch_an_counts <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Corrected/MBatch_AN/Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
mbatch_mp_counts <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Corrected/MBatch_MP/Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
uncorrected_counts <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Corrected/Uncorrected/Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)

dataset_batchqc <- matrix(c(rep("GLDS-47", 6), rep("GLDS-48", 14), rep("GLDS-137", 12), rep("GLDS-168", 28), rep("GLDS-173", 4), rep("GLDS-242", 9),rep("GLDS-245", 39)))


# batchQC(uncorrected_counts, batch=dataset_batchqc, report_file = "uncorrected_liver_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")
# batchQC(combat_counts, batch=dataset_batchqc, report_file = "combat_liver_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")
# batchQC(combat_seq_counts, batch=dataset_batchqc, report_file = "combat_seq_liver_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")
# batchQC(deseq2_counts, batch=dataset_batchqc, report_file = "deseq2_liver_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")
# batchQC(mbatch_eb_counts, batch=dataset_batchqc, report_file = "mbatch_eb_liver_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")
# batchQC(mbatch_an_counts, batch=dataset_batchqc, report_file = "mbatch_an_liver_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")
batchQC(mbatch_mp_counts, batch=dataset_batchqc, report_file = "mbatch_mp_liver_report.html", report_dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BatchQC/reports")


