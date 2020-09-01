# Generates corrected counts with various methods
# Covariates were not used for corrections
# Finsam Samson

library(devtools)
library(Biobase)
library(sva)
library(tximport)
library(DESeq2)
library(tidyverse)
library(ggfortify)
library(ggplot2)
library(MBatch)
library(limma)
library(parallel)
library(BiocParallel)



# Import unnormalized raw counts as data frames
# Extract flight and ground control samples, combine
raw_47 <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-47/GLDS-47_rna_seq_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
raw_47 <- raw_47[(!grepl("ERCC", rownames(raw_47))),]
raw_47_flt <- as.matrix(raw_47[,grepl("FLT", names(raw_47))])
raw_47_gc <- as.matrix(raw_47[,grepl("GC", names(raw_47))])
raw_47_flt_gc <- cbind(raw_47_flt, raw_47_gc)

raw_48 <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-48/GLDS-48_rna_seq_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
raw_48 <- raw_48[(!grepl("ERCC", rownames(raw_48))),]
raw_48_flt <- as.matrix(raw_48[,grepl("FLT", names(raw_48))])
raw_48_gc <- as.matrix(raw_48[,grepl("GC", names(raw_48))])
raw_48_flt_gc <- cbind(raw_48_flt, raw_48_gc)

raw_137 <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-137/GLDS-137_rna_seq_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
raw_137 <- raw_137[(!grepl("ERCC", rownames(raw_137))),]
raw_137_flt <- as.matrix(raw_137[,grepl("FLT", names(raw_137))])
raw_137_gc <- as.matrix(raw_137[,grepl("GC", names(raw_137))])
raw_137_flt_gc <- cbind(raw_137_flt, raw_137_gc)

raw_168 <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-168/GLDS-168_rna_seq_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
raw_168 <- raw_168[(!grepl("ERCC", rownames(raw_168))),]
raw_168_flt <- as.matrix(raw_168[,grepl("FLT", names(raw_168))])
raw_168_gc <- as.matrix(raw_168[,grepl("GC", names(raw_168))])
raw_168_flt_gc <- cbind(raw_168_flt, raw_168_gc)

## new sample table, metadata table for GLDS-168
sample_table_168 <- matrix(c(rep("FLT", 14), rep("GC", 14)))
rownames(sample_table_168) <- colnames(raw_168_flt_gc)
metadata_168 <- sample_table_168
colnames(metadata_168) <- "condition"
colnames(sample_table_168) <- "Group"

write.csv(raw_168_flt_gc,file.path("~/Documents/SLSTP/Raw/GLDS-168/GLDS-168_input_files", "Unnormalized_Counts.csv"), row.names = TRUE)
write.csv(sample_table_168,file.path("~/Documents/SLSTP/Raw/GLDS-168/GLDS-168_input_files", "GLDS-168-FLT-GC_sample_group.csv"), row.names = TRUE)
write.csv(metadata_168,file.path("~/Documents/SLSTP/Raw/GLDS-168/GLDS-168_input_files", "GLDS-168-FLT-GC_metadata.csv"), row.names = TRUE)

raw_173 <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-173/GLDS-173_rna_seq_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
raw_173 <- raw_173[(!grepl("ERCC", rownames(raw_173))),]
raw_173_flt <- as.matrix(raw_173[,grepl("FLT", names(raw_173))])
raw_173_gc <- as.matrix(raw_173[,grepl("GC", names(raw_173))])
raw_173_flt_gc <- cbind(raw_173_flt, raw_173_gc)

raw_242 <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-242/GLDS-242_rna_seq_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
raw_242 <- raw_242[(!grepl("ERCC", rownames(raw_242))),]
raw_242_flt <- as.matrix(raw_242[,grepl("FLT", names(raw_242))])
raw_242_gc <- as.matrix(raw_242[,grepl("GC", names(raw_242))])
raw_242_flt_gc <- cbind(raw_242_flt, raw_242_gc)

raw_245 <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Raw/GLDS-245/GLDS-245_rna_seq_Unnormalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
raw_245 <- raw_245[(!grepl("ERCC", rownames(raw_245))),]
raw_245_flt <- as.matrix(raw_245[,grepl("FLT", names(raw_245))])
raw_245_gc <- as.matrix(raw_245[,grepl("GC", names(raw_245))])
raw_245_flt_gc <- cbind(raw_245_flt, raw_245_gc)

# Merge all datasets
raw_liver <- cbind(raw_47_flt_gc, raw_48_flt_gc, raw_137_flt_gc, raw_168_flt_gc, raw_173_flt_gc, raw_242_flt_gc, raw_245_flt_gc)

# Create dataset/sample/condition tables for corrections
dataset_int <- c(rep.int(47, 6), rep.int(48, 14), rep.int(137, 12), rep.int(168, 28), rep.int(173, 4), rep.int(242, 9),rep.int(245, 39))
dataset_str <- data.frame(c(rep("GLDS-47", 6), rep("GLDS-48", 14), rep("GLDS-137", 12), rep("GLDS-168", 28), rep("GLDS-173", 4), rep("GLDS-242", 9),rep("GLDS-245", 39)))
colnames(dataset_str) <- "Group"
# 1 = FLT, 2 = GC
condition_int <- c(rep.int(1, 3), rep.int(2, 3), rep.int(1, 7), rep.int(2, 7), rep.int(1, 6), rep.int(2, 6), rep.int(1, 14), rep.int(2, 14), rep.int(1, 2), rep.int(2, 2), rep.int(1, 5), rep.int(2, 4), rep.int(1, 20), rep.int(2, 19))
condition_str <- c(rep("FLT", 3), rep("GC", 3), rep("FLT", 7), rep("GC", 7), rep("FLT", 6), rep("GC", 6), rep("FLT", 14), rep("GC", 14), rep("FLT", 2), rep("GC", 2), rep("FLT", 5), rep("GC", 4), rep("FLT", 20), rep("GC", 19))


### Apply correction methods

## ComBat
combat_liver <- ComBat(dat=raw_liver, batch=dataset_int, mod=NULL, par.prior=TRUE, prior.plots=FALSE)



## ComBat_seq
# only worked without covariates
combat_seq_liver <- ComBat_seq(counts=raw_liver, batch=dataset_int)



## MBatch - EB
# Define paths to input/output directories
work_dir="~/Documents/SLSTP/Raw/MBatch_Input"
# The out_dir must be empty - this is where the corrected output table will be written
out_dir="~/Documents/SLSTP/Processed/MBatch_Output"

# Pull in data and metadata files
# Raw counts table, sample_names as column headers and row names are geneIDs
theData <- (file.path(work_dir,"all_liver_counts.tsv"))
# Table with the sample_names as row names and column headers are the batch_types, each sample has a defined batch for each batch_type
theBatches <- (file.path(work_dir,"theBatches.tsv"))
# Define which batch_type you want to correct for
theBatchType="GLDS"
# Table with the sample_names as row names and column headers are the covariate_types
theCovariates <- (file.path(work_dir,"theCovariates.tsv"))

# Load the data in MBatch format (used in both EB and AN)
myData <- mbatchLoadFiles(theData, theBatches)

# nonparametric priors did not finish in ~6 hours
eb_liver <- EB_withParametricPriors(theBeaData=myData,
                                         theBatchIdsNotToCorrect=c(""),
                                         theDoCheckPlotsFlag=FALSE,
                                         theBatchType=theBatchType,
                                         theThreads=1,
                                         thePath=out_dir,
                                         theWriteToFile=FALSE)


## MBatch - AN
an_liver <- AN_Adjusted(theBeaData=myData,
                             theBatchType=theBatchType,
                             thePath=out_dir,
                             theWriteToFile=FALSE)



## MBatch - MP
output_dir <- "~/Documents/SLSTP/Processed/MBatch_Output/mp_file.tsv"
empty <- data.frame()
# Create MBatch BEA Object for MP
bea_liver <- new("BEA_DATA", raw_liver, dataset_str, empty)
mp_liver <- MP_Overall(theBeaData=bea_liver, thePath=output_dir, theWriteToFile=FALSE)



## Prepare all tables for export
# for dge export
combat_dge <- combat_liver
combat_seq_dge <- combat_seq_liver
deseq2_dge <- raw_liver
mbatch_eb_dge <- eb_liver
mbatch_an_dge <- an_liver
mbatch_mp_dge <- mp_liver
uncorrected <- raw_liver

# Convert negatives to 0
combat_dge[combat_dge<0] <- 0
combat_seq_dge[combat_seq_dge<0] <- 0
deseq2_dge[deseq2_dge<0] <- 0
mbatch_eb_dge[mbatch_eb_dge<0] <- 0
mbatch_an_dge[mbatch_an_dge<0] <- 0
mbatch_an_dge[is.nan(mbatch_an_dge)] <- 0
mbatch_mp_dge[mbatch_mp_dge<0] <- 0
mbatch_mp_dge[is.nan(mbatch_mp_dge)] <- 0
mbatch_mp_dge <- mbatch_mp_dge + 1
uncorrected[uncorrected<0] <- 0


## Export sample groups/metadata
# For DESeq2 DGE
export_groups <- matrix(c(rep("47_FLT", 3), rep("47_GC", 3), rep("48_FLT", 7), rep("48_GC", 7), rep("137_FLT", 6), rep("137_GC", 6), rep("168_FLT", 14), rep("168_GC", 14), rep("173_FLT", 2), rep("173_GC", 2), rep("242_FLT", 5), rep("242_GC", 4), rep("245_FLT", 20), rep("245_GC", 19)))
rownames(export_groups) = colnames(raw_liver)
colnames(export_groups) = 'Group'
# write.csv(export_groups,file.path("~/Documents/SLSTP/Corrected/Uncorrected", "all_liver_sample_group.csv"), row.names = TRUE)

condition_exp <- matrix(c(rep("FLT", 3), rep("GC", 3), rep("FLT", 7), rep("GC", 7), rep("FLT", 6), rep("GC", 6), rep("FLT", 14), rep("GC", 14), rep("FLT", 2), rep("GC", 2), rep("FLT", 5), rep("GC", 4), rep("FLT", 20), rep("GC", 19)))
colnames(condition_exp) = 'condition'
dataset_exp <- matrix(c(rep("GLDS-47", 6), rep("GLDS-48", 14), rep("GLDS-137", 12), rep("GLDS-168", 28), rep("GLDS-173", 4), rep("GLDS-242", 9),rep("GLDS-245", 39)))
colnames(dataset_exp) = 'dataset'
export_meta <- cbind(condition_exp, dataset_exp)
rownames(export_meta) = colnames(raw_liver)
# write.csv(export_meta,file.path("~/Documents/SLSTP/Corrected/Uncorrected", "all_liver_metadata.csv"), row.names = TRUE)


## Alternate groupings for DESeq2 normalization
# condition_exp <- matrix(c(rep("FLT", 3), rep("GC", 3), rep("FLT", 7), rep("GC", 7), rep("FLT", 6), rep("GC", 6), rep("FLT", 14), rep("GC", 14), rep("FLT", 2), rep("GC", 2), rep("FLT", 5), rep("GC", 4), rep("FLT", 20), rep("GC", 19)))
condition_exp <- matrix(c(rep("47_FLT", 3), rep("47_GC", 3), rep("48_FLT", 7), rep("48_GC", 7), rep("137_FLT", 6), rep("137_GC", 6), rep("168_FLT", 14), rep("168_GC", 14), rep("173_FLT", 2), rep("173_GC", 2), rep("242_FLT", 5), rep("242_GC", 4), rep("245_FLT", 20), rep("245_GC", 19)))
colnames(condition_exp) = 'condition'
rownames(condition_exp) = colnames(raw_liver)
# write.csv(condition_exp,file.path("~/Documents/SLSTP/Corrected/DESeq2", "all_liver_metadata.csv"), row.names = TRUE)



## EXPORT counts
# write.csv(combat_dge,file.path("~/Documents/SLSTP/Corrected/ComBat", "Unnormalized_Counts.csv"), row.names = TRUE)
# write.csv(combat_seq_dge,file.path("~/Documents/SLSTP/Corrected/ComBat_seq", "Unnormalized_Counts.csv"), row.names = TRUE)
# write.csv(deseq2_dge,file.path("~/Documents/SLSTP/Corrected/DESeq2", "Unnormalized_Counts.csv"), row.names = TRUE)
# write.csv(mbatch_eb_dge,file.path("~/Documents/SLSTP/Corrected/MBatch_EB", "Unnormalized_Counts.csv"), row.names = TRUE)
# write.csv(mbatch_an_dge,file.path("~/Documents/SLSTP/Corrected/MBatch_AN", "Unnormalized_Counts.csv"), row.names = TRUE)
# write.csv(mbatch_mp_dge,file.path("~/Documents/SLSTP/Corrected/MBatch_MP", "Unnormalized_Counts.csv"), row.names = TRUE)
# write.csv(uncorrected,file.path("~/Documents/SLSTP/Corrected/Uncorrected", "Unnormalized_Counts.csv"), row.names = TRUE)

