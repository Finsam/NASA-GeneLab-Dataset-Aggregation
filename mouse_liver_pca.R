# run PCAs on liver data
# Finsam Samson - 081520

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

## import uncorrected/corrected counts tables
combat_counts <- as.matrix(read.csv(Sys.glob(file.path("~/Documents/SLSTP/Corrected/ComBat/Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))
combat_seq_counts <- as.matrix(read.csv(Sys.glob(file.path("~/Documents/SLSTP/Corrected/ComBat_seq/Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))
deseq2_counts <- as.matrix(read.csv(Sys.glob(file.path("~/Documents/SLSTP/Processed/All_Liver/DESeq2/DSeq2_NormCounts/Normalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))
mbatch_eb_counts <- as.matrix(read.csv(Sys.glob(file.path("~/Documents/SLSTP/Corrected/MBatch_EB/Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))
mbatch_an_counts <- as.matrix(read.csv(Sys.glob(file.path("~/Documents/SLSTP/Corrected/MBatch_AN/Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))
mbatch_mp_counts <- as.matrix(read.csv(Sys.glob(file.path("~/Documents/SLSTP/Corrected/MBatch_MP/Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))
uncorrected_counts <- as.matrix(read.csv(Sys.glob(file.path("~/Documents/SLSTP/Corrected/Uncorrected/Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))


# PCA Generation
# Get log2 counts
uncorrected_in <- log2(uncorrected_counts) 
combat_in <- log2(combat_counts) 
combat_seq_in <- log2(combat_seq_counts) 
deseq2_in <- log2(deseq2_counts) 
mbatch_eb_in <- log2(mbatch_eb_counts) 
mbatch_an_in <- log2(mbatch_an_counts) 
mbatch_mp_in <- log2(mbatch_mp_counts) 

# Replace all NAN values with 0
uncorrected_in[is.nan(uncorrected_in)] <- 0
combat_in[is.nan(combat_in)] <- 0 
combat_seq_in[is.nan(combat_seq_in)] <- 0
deseq2_in[is.nan(deseq2_in)] <- 0
mbatch_eb_in[is.nan(mbatch_eb_in)] <- 0
mbatch_an_in[is.nan(mbatch_an_in)] <- 0
mbatch_mp_in[is.nan(mbatch_mp_in)] <- 0

uncorrected_in[is.infinite(uncorrected_in)] <- 0
combat_in[is.infinite(combat_in)] <- 0 
combat_seq_in[is.infinite(combat_seq_in)] <- 0
deseq2_in[is.infinite(deseq2_in)] <- 0
mbatch_eb_in[is.infinite(mbatch_eb_in)] <- 0
mbatch_an_in[is.infinite(mbatch_an_in)] <- 0
mbatch_mp_in[is.infinite(mbatch_mp_in)] <- 0

## Run PCA
uncorrected_pca <- prcomp(t(uncorrected_in), scale = FALSE)
combat_pca <- prcomp(t(combat_in), scale = FALSE)
combat_seq_pca <- prcomp(t(combat_seq_in), scale = FALSE)
deseq2_pca <- prcomp(t(deseq2_in), scale = FALSE)
mbatch_eb_pca <- prcomp(t(mbatch_eb_in), scale = FALSE)
mbatch_an_pca <- prcomp(t(mbatch_an_in), scale = FALSE)
mbatch_mp_pca <- prcomp(t(mbatch_mp_in), scale = FALSE)


## import bigger metadata file
metadata_no_d <- read_tsv(Sys.glob(file.path("~/Downloads/theCovariates.txt")), col_names = TRUE)

# Create PCA groups from metadata file
# euthanasia_pca <- metadata[,"EuthanasiaLoc"]
# spaceflight_pca <- metadata[,"Spaceflight"]
# rna_prep_pca <- metadata[,"RNAprep"]
# sex_pca <- metadata[,"Sex"]
# strain_pca <- metadata[,"Strain"]
dataset_pca <- data.frame(c(rep("GLDS-47", 6), rep("GLDS-48", 14), rep("GLDS-137", 12), rep("GLDS-168", 28), rep("GLDS-173", 4), rep("GLDS-242", 9),rep("GLDS-245", 39)))
colnames(dataset_pca) <- "Dataset"
metadata_pca <- cbind(metadata_no_d, dataset_pca)


## Plot PCAs
size_var <- 3
alpha_var <- 0.6

autoplot(uncorrected_pca, data=metadata_pca, colour='Dataset', shape='Spaceflight', size=size_var, alpha=alpha_var) + scale_shape_manual(values=c(17,16)) + theme_classic()
autoplot(combat_pca, data=metadata_pca, colour='Dataset', shape='Spaceflight', size=size_var, alpha=alpha_var) + scale_shape_manual(values=c(17,16)) + theme_classic()
autoplot(combat_seq_pca, data=metadata_pca, colour='Dataset', shape='Spaceflight', size=size_var, alpha=alpha_var) + scale_shape_manual(values=c(17,16)) + theme_classic()
autoplot(deseq2_pca, data=metadata_pca, colour='Dataset', shape='Spaceflight', size=size_var, alpha=alpha_var) + scale_shape_manual(values=c(17,16)) + theme_classic()
autoplot(mbatch_eb_pca, data=metadata_pca, colour='Dataset', shape='Spaceflight', size=size_var, alpha=alpha_var) + scale_shape_manual(values=c(17,16)) + theme_classic()
autoplot(mbatch_an_pca, data=metadata_pca, colour='Dataset', shape='Spaceflight', size=size_var, alpha=alpha_var) + scale_shape_manual(values=c(17,16)) + theme_classic()
autoplot(mbatch_mp_pca, data=metadata_pca, colour='Dataset', shape='Spaceflight', size=size_var, alpha=alpha_var) + scale_shape_manual(values=c(17,16)) + theme_classic()
