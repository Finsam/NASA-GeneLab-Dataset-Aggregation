# script to filter DGE tables
# Finsam Samson
#
#
# After running:
## total_table contains count of total DEGs
## match_table contains count of DEGs that match those of original dataset
## extra_table contains count of DEGs that don't match those of original

library(BatchQC)
library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)
library(tximport)
library(DESeq2)
library(tidyverse)
library(ggfortify)
library(ggplot2)
library(MBatch)
library(limma)
library(edgeR)
library(parallel)
library(BiocParallel)
library(dplyr)

## Import the DGE tables of interest
## In this case I wanted tables with FLT v. GC groups
# dge_47 <- read.csv(Sys.glob(file.path("~/Downloads/GLDS-47_rna_seq_differential_expression.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# dge_48 <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Processed/GLDS-48/Uncorrected/DSeq2_DGE/differential_expression.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# dge_137 <- read.csv(Sys.glob(file.path("~/Downloads/GLDS-137_rna_seq_differential_expression.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# dge_168 <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Processed/GLDS-168/Uncorrected/DSeq2_DGE/differential_expression.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# dge_173 <- read.csv(Sys.glob(file.path("~/Downloads/GLDS-173_rna_seq_differential_expression.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# dge_242 <- read.csv(Sys.glob(file.path("~/Downloads/GLDS-242_rna_seq_differential_expression.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# dge_245 <- read.csv(Sys.glob(file.path("~/Downloads/GLDS-245_differential_expression_ISST_LAR_grouped.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)


## GLDS unfiltered subsets
trim_47 <- select(dge_47, "Log2fc_.FLT.v.GC.", "P.value_.FLT.v.GC.", "Adj.p.value_.FLT.v.GC.")
trim_48 <- select(dge_48, "Log2fc_.FLT.v.GC.", "P.value_.FLT.v.GC.", "Adj.p.value_.FLT.v.GC.")
# trim_48_C <- select(dge_48, "Log2fc_.FLT_C.v.GC_C.", "P.value_.FLT_C.v.GC_C.", "Adj.p.value_.FLT_C.v.GC_C.")
# trim_48_I <- select(dge_48, "Log2fc_.FLT_I.v.GC_I.", "P.value_.FLT_I.v.GC_I.", "Adj.p.value_.FLT_I.v.GC_I.")
trim_137 <- select(dge_137, "Log2fc_.FLT.v.GC.", "P.value_.FLT.v.GC.", "Adj.p.value_.FLT.v.GC.")
trim_168 <- select(dge_168, "Log2fc_.FLT.v.GC.", "P.value_.FLT.v.GC.", "Adj.p.value_.FLT.v.GC.")
# With/without ERCC, or pick one; doing separate rodent missions separately
# NEED NEW DGE: trim_168 <- select(dge_168, "Log2fc_.FLT.v.GC.", "P.value_.FLT.v.GC.", "Adj.p.value_.FLT.v.GC.")
# trim_168_1 <- select(dge_168, "Log2fc_.RR1_FLT_noERCC.v.RR1_GC_noERCC.", "P.value_.RR1_FLT_noERCC.v.RR1_GC_noERCC.", "Adj.p.value_.RR1_FLT_noERCC.v.RR1_GC_noERCC.")
# trim_168_3 <- select(dge_168, "Log2fc_.RR3_FLT_wERCC.v.RR3_GC_wERCC.", "P.value_.RR3_FLT_wERCC.v.RR3_GC_wERCC.", "Adj.p.value_.RR3_FLT_wERCC.v.RR3_GC_wERCC.")
trim_173 <- select(dge_173, "Log2fc_.FLT.v.GC.", "P.value_.FLT.v.GC.", "Adj.p.value_.FLT.v.GC.")
trim_242 <- select(dge_242, "Log2fc_.FLT_C1.v.GC_C2.", "P.value_.FLT_C1.v.GC_C2.", "Adj.p.value_.FLT_C1.v.GC_C2.")
trim_245 <- select(dge_245, "Log2fc_.FLT.v.GC.", "P.value_.FLT.v.GC.", "Adj.p.value_.FLT.v.GC.")

## GLDS filtered by p
p_47 <- trim_47[trim_47$Adj.p.value_.FLT.v.GC.<0.05,]
p_47 <- p_47[is.finite(p_47$Adj.p.value_.FLT.v.GC.),]
p_48 <- trim_48[trim_48$Adj.p.value_.FLT.v.GC.<0.05,]
p_48 <- p_48[is.finite(p_48$Adj.p.value_.FLT.v.GC.),]
# p_48_C <- trim_48_C[trim_48_C$Adj.p.value_.FLT_C.v.GC_C.<0.05,]
# p_48_C <- p_48_C[is.finite(p_48_C$Adj.p.value_.FLT_C.v.GC_C.),]
# p_48_I <- trim_48_I[trim_48_I$Adj.p.value_.FLT_I.v.GC_I.<0.05,]
# p_48_I <- p_48_I[is.finite(p_48_I$Adj.p.value_.FLT_I.v.GC_I.),]
p_137 <- trim_137[trim_137$Adj.p.value_.FLT.v.GC.<0.05,]
p_137 <- p_137[is.finite(p_137$Adj.p.value_.FLT.v.GC.),]
p_168 <- trim_168[trim_168$Adj.p.value_.FLT.v.GC.<0.05,]
p_168 <- p_168[is.finite(p_168$Adj.p.value_.FLT.v.GC.),]
# p_168_1 <- trim_168_1[trim_168_1$Adj.p.value_.RR1_FLT_noERCC.v.RR1_GC_noERCC.<0.05,]
# p_168_1 <- p_168_1[is.finite(p_168_1$Adj.p.value_.RR1_FLT_noERCC.v.RR1_GC_noERCC.),]
# p_168_3 <- trim_168_3[trim_168_3$Adj.p.value_.RR3_FLT_wERCC.v.RR3_GC_wERCC.<0.05,]
# p_168_3 <- p_168_3[is.finite(p_168_3$Adj.p.value_.RR3_FLT_wERCC.v.RR3_GC_wERCC.),]
p_173 <- trim_173[trim_173$Adj.p.value_.FLT.v.GC.<0.05,]
p_173 <- p_173[is.finite(p_173$Adj.p.value_.FLT.v.GC.),]
p_242 <- trim_242[trim_242$Adj.p.value_.FLT_C1.v.GC_C2.<0.05,]
p_242 <- p_242[is.finite(p_242$Adj.p.value_.FLT_C1.v.GC_C2.),]
p_245 <- trim_245[trim_245$Adj.p.value_.FLT.v.GC.<0.05,]
p_245 <- p_245[is.finite(p_245$Adj.p.value_.FLT.v.GC.),]

# For filtering:
# Log2fc_(FLT & GLDS-173)v(GC & GLDS-173)
# P.value_(FLT & GLDS-173)v(GC & GLDS-173)
# Adj.p.value_(FLT & GLDS-173)v(GC & GLDS-173)
# Log2fc_.GC...GLDS.47.v.GC...GLDS.242.

# Import corrected dge tables
# uncorrected_dge <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Processed/All_Liver/Uncorrected/DSeq2_DGE/differential_expression.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# combat_dge <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Processed/All_Liver/ComBat/DSeq2_DGE/differential_expression.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# combat_seq_dge <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Processed/All_Liver/ComBat_seq/DSeq2_DGE/differential_expression.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# deseq2_dge <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Processed/All_Liver/DESeq2/DSeq2_DGE/differential_expression.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# mbatch_eb_dge <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Processed/All_Liver/MBatch_EB/DSeq2_DGE/differential_expression.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# mbatch_an_dge <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Processed/All_Liver/MBatch_AN/DSeq2_DGE/differential_expression.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
# mbatch_mp_dge <- read.csv(Sys.glob(file.path("~/Documents/SLSTP/Processed/All_Liver/MBatch_MP/DSeq2_DGE/differential_expression.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)



## TRIM all columns except for FLT v. GC columns

## uncorrected unfiltered subsets
uncorrected_trim_47 <- select(uncorrected_dge, "Log2fc_.FLT...GLDS.47.v.GC...GLDS.47.", "P.value_.FLT...GLDS.47.v.GC...GLDS.47.", "Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.")
uncorrected_trim_48 <- select(uncorrected_dge, "Log2fc_.FLT...GLDS.48.v.GC...GLDS.48.", "P.value_.FLT...GLDS.48.v.GC...GLDS.48.", "Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.")
uncorrected_trim_137 <- select(uncorrected_dge, "Log2fc_.FLT...GLDS.137.v.GC...GLDS.137.", "P.value_.FLT...GLDS.137.v.GC...GLDS.137.", "Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.")
uncorrected_trim_168 <- select(uncorrected_dge, "Log2fc_.FLT...GLDS.168.v.GC...GLDS.168.", "P.value_.FLT...GLDS.168.v.GC...GLDS.168.", "Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.")
uncorrected_trim_173 <- select(uncorrected_dge, "Log2fc_.FLT...GLDS.173.v.GC...GLDS.173.", "P.value_.FLT...GLDS.173.v.GC...GLDS.173.", "Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.")
uncorrected_trim_242 <- select(uncorrected_dge, "Log2fc_.FLT...GLDS.242.v.GC...GLDS.242.", "P.value_.FLT...GLDS.242.v.GC...GLDS.242.", "Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.")
uncorrected_trim_245 <- select(uncorrected_dge, "Log2fc_.FLT...GLDS.245.v.GC...GLDS.245.", "P.value_.FLT...GLDS.245.v.GC...GLDS.245.", "Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.")

## ComBat unfiltered subsets
combat_trim_47 <- select(combat_dge, "Log2fc_.FLT...GLDS.47.v.GC...GLDS.47.", "P.value_.FLT...GLDS.47.v.GC...GLDS.47.", "Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.")
combat_trim_48 <- select(combat_dge, "Log2fc_.FLT...GLDS.48.v.GC...GLDS.48.", "P.value_.FLT...GLDS.48.v.GC...GLDS.48.", "Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.")
combat_trim_137 <- select(combat_dge, "Log2fc_.FLT...GLDS.137.v.GC...GLDS.137.", "P.value_.FLT...GLDS.137.v.GC...GLDS.137.", "Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.")
combat_trim_168 <- select(combat_dge, "Log2fc_.FLT...GLDS.168.v.GC...GLDS.168.", "P.value_.FLT...GLDS.168.v.GC...GLDS.168.", "Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.")
combat_trim_173 <- select(combat_dge, "Log2fc_.FLT...GLDS.173.v.GC...GLDS.173.", "P.value_.FLT...GLDS.173.v.GC...GLDS.173.", "Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.")
combat_trim_242 <- select(combat_dge, "Log2fc_.FLT...GLDS.242.v.GC...GLDS.242.", "P.value_.FLT...GLDS.242.v.GC...GLDS.242.", "Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.")
combat_trim_245 <- select(combat_dge, "Log2fc_.FLT...GLDS.245.v.GC...GLDS.245.", "P.value_.FLT...GLDS.245.v.GC...GLDS.245.", "Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.")

## ComBat_seq unfiltered subsets
combat_seq_trim_47 <- select(combat_seq_dge, "Log2fc_.FLT...GLDS.47.v.GC...GLDS.47.", "P.value_.FLT...GLDS.47.v.GC...GLDS.47.", "Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.")
combat_seq_trim_48 <- select(combat_seq_dge, "Log2fc_.FLT...GLDS.48.v.GC...GLDS.48.", "P.value_.FLT...GLDS.48.v.GC...GLDS.48.", "Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.")
combat_seq_trim_137 <- select(combat_seq_dge, "Log2fc_.FLT...GLDS.137.v.GC...GLDS.137.", "P.value_.FLT...GLDS.137.v.GC...GLDS.137.", "Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.")
combat_seq_trim_168 <- select(combat_seq_dge, "Log2fc_.FLT...GLDS.168.v.GC...GLDS.168.", "P.value_.FLT...GLDS.168.v.GC...GLDS.168.", "Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.")
combat_seq_trim_173 <- select(combat_seq_dge, "Log2fc_.FLT...GLDS.173.v.GC...GLDS.173.", "P.value_.FLT...GLDS.173.v.GC...GLDS.173.", "Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.")
combat_seq_trim_242 <- select(combat_seq_dge, "Log2fc_.FLT...GLDS.242.v.GC...GLDS.242.", "P.value_.FLT...GLDS.242.v.GC...GLDS.242.", "Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.")
combat_seq_trim_245 <- select(combat_seq_dge, "Log2fc_.FLT...GLDS.245.v.GC...GLDS.245.", "P.value_.FLT...GLDS.245.v.GC...GLDS.245.", "Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.")

##
## Waiting to do DESeq2
##

## MBatch_EB unfiltered subsets
mbatch_eb_trim_47 <- select(mbatch_eb_dge, "Log2fc_.FLT...GLDS.47.v.GC...GLDS.47.", "P.value_.FLT...GLDS.47.v.GC...GLDS.47.", "Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.")
mbatch_eb_trim_48 <- select(mbatch_eb_dge, "Log2fc_.FLT...GLDS.48.v.GC...GLDS.48.", "P.value_.FLT...GLDS.48.v.GC...GLDS.48.", "Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.")
mbatch_eb_trim_137 <- select(mbatch_eb_dge, "Log2fc_.FLT...GLDS.137.v.GC...GLDS.137.", "P.value_.FLT...GLDS.137.v.GC...GLDS.137.", "Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.")
mbatch_eb_trim_168 <- select(mbatch_eb_dge, "Log2fc_.FLT...GLDS.168.v.GC...GLDS.168.", "P.value_.FLT...GLDS.168.v.GC...GLDS.168.", "Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.")
mbatch_eb_trim_173 <- select(mbatch_eb_dge, "Log2fc_.FLT...GLDS.173.v.GC...GLDS.173.", "P.value_.FLT...GLDS.173.v.GC...GLDS.173.", "Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.")
mbatch_eb_trim_242 <- select(mbatch_eb_dge, "Log2fc_.FLT...GLDS.242.v.GC...GLDS.242.", "P.value_.FLT...GLDS.242.v.GC...GLDS.242.", "Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.")
mbatch_eb_trim_245 <- select(mbatch_eb_dge, "Log2fc_.FLT...GLDS.245.v.GC...GLDS.245.", "P.value_.FLT...GLDS.245.v.GC...GLDS.245.", "Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.")

## MBatch_AN unfiltered subsets
mbatch_an_trim_47 <- select(mbatch_an_dge, "Log2fc_.FLT...GLDS.47.v.GC...GLDS.47.", "P.value_.FLT...GLDS.47.v.GC...GLDS.47.", "Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.")
mbatch_an_trim_48 <- select(mbatch_an_dge, "Log2fc_.FLT...GLDS.48.v.GC...GLDS.48.", "P.value_.FLT...GLDS.48.v.GC...GLDS.48.", "Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.")
mbatch_an_trim_137 <- select(mbatch_an_dge, "Log2fc_.FLT...GLDS.137.v.GC...GLDS.137.", "P.value_.FLT...GLDS.137.v.GC...GLDS.137.", "Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.")
mbatch_an_trim_168 <- select(mbatch_an_dge, "Log2fc_.FLT...GLDS.168.v.GC...GLDS.168.", "P.value_.FLT...GLDS.168.v.GC...GLDS.168.", "Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.")
mbatch_an_trim_173 <- select(mbatch_an_dge, "Log2fc_.FLT...GLDS.173.v.GC...GLDS.173.", "P.value_.FLT...GLDS.173.v.GC...GLDS.173.", "Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.")
mbatch_an_trim_242 <- select(mbatch_an_dge, "Log2fc_.FLT...GLDS.242.v.GC...GLDS.242.", "P.value_.FLT...GLDS.242.v.GC...GLDS.242.", "Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.")
mbatch_an_trim_245 <- select(mbatch_an_dge, "Log2fc_.FLT...GLDS.245.v.GC...GLDS.245.", "P.value_.FLT...GLDS.245.v.GC...GLDS.245.", "Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.")

## MBatch_MP unfiltered subsets
mbatch_mp_trim_47 <- select(mbatch_mp_dge, "Log2fc_.FLT...GLDS.47.v.GC...GLDS.47.", "P.value_.FLT...GLDS.47.v.GC...GLDS.47.", "Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.")
mbatch_mp_trim_48 <- select(mbatch_mp_dge, "Log2fc_.FLT...GLDS.48.v.GC...GLDS.48.", "P.value_.FLT...GLDS.48.v.GC...GLDS.48.", "Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.")
mbatch_mp_trim_137 <- select(mbatch_mp_dge, "Log2fc_.FLT...GLDS.137.v.GC...GLDS.137.", "P.value_.FLT...GLDS.137.v.GC...GLDS.137.", "Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.")
mbatch_mp_trim_168 <- select(mbatch_mp_dge, "Log2fc_.FLT...GLDS.168.v.GC...GLDS.168.", "P.value_.FLT...GLDS.168.v.GC...GLDS.168.", "Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.")
mbatch_mp_trim_173 <- select(mbatch_mp_dge, "Log2fc_.FLT...GLDS.173.v.GC...GLDS.173.", "P.value_.FLT...GLDS.173.v.GC...GLDS.173.", "Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.")
mbatch_mp_trim_242 <- select(mbatch_mp_dge, "Log2fc_.FLT...GLDS.242.v.GC...GLDS.242.", "P.value_.FLT...GLDS.242.v.GC...GLDS.242.", "Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.")
mbatch_mp_trim_245 <- select(mbatch_mp_dge, "Log2fc_.FLT...GLDS.245.v.GC...GLDS.245.", "P.value_.FLT...GLDS.245.v.GC...GLDS.245.", "Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.")



## Filter to include genes with adj. p < 0.05
# uncorrected 
uncorrected_p_47 <- uncorrected_trim_47[uncorrected_trim_47$Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.<0.05,]
uncorrected_p_47 <- uncorrected_p_47[is.finite(uncorrected_p_47$Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.),]
uncorrected_p_48 <- uncorrected_trim_48[uncorrected_trim_48$Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.<0.05,]
uncorrected_p_48 <- uncorrected_p_48[is.finite(uncorrected_p_48$Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.),]
uncorrected_p_137 <- uncorrected_trim_137[uncorrected_trim_137$Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.<0.05,]
uncorrected_p_137 <- uncorrected_p_137[is.finite(uncorrected_p_137$Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.),]
uncorrected_p_168 <- uncorrected_trim_168[uncorrected_trim_168$Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.<0.05,]
uncorrected_p_168 <- uncorrected_p_168[is.finite(uncorrected_p_168$Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.),]
uncorrected_p_173 <- uncorrected_trim_173[uncorrected_trim_173$Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.<0.05,]
uncorrected_p_173 <- uncorrected_p_173[is.finite(uncorrected_p_173$Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.),]
uncorrected_p_242 <- uncorrected_trim_242[uncorrected_trim_242$Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.<0.05,]
uncorrected_p_242 <- uncorrected_p_242[is.finite(uncorrected_p_242$Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.),]
uncorrected_p_245 <- uncorrected_trim_245[uncorrected_trim_245$Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.<0.05,]
uncorrected_p_245 <- uncorrected_p_245[is.finite(uncorrected_p_245$Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.),]

# combat
combat_p_47 <- combat_trim_47[combat_trim_47$Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.<0.05,]
combat_p_47 <- combat_p_47[is.finite(combat_p_47$Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.),]
combat_p_48 <- combat_trim_48[combat_trim_48$Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.<0.05,]
combat_p_48 <- combat_p_48[is.finite(combat_p_48$Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.),]
combat_p_137 <- combat_trim_137[combat_trim_137$Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.<0.05,]
combat_p_137 <- combat_p_137[is.finite(combat_p_137$Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.),]
combat_p_168 <- combat_trim_168[combat_trim_168$Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.<0.05,]
combat_p_168 <- combat_p_168[is.finite(combat_p_168$Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.),]
combat_p_173 <- combat_trim_173[combat_trim_173$Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.<0.05,]
combat_p_173 <- combat_p_173[is.finite(combat_p_173$Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.),]
combat_p_242 <- combat_trim_242[combat_trim_242$Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.<0.05,]
combat_p_242 <- combat_p_242[is.finite(combat_p_242$Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.),]
combat_p_245 <- combat_trim_245[combat_trim_245$Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.<0.05,]
combat_p_245 <- combat_p_245[is.finite(combat_p_245$Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.),]

# combat_seq
combat_seq_p_47 <- combat_seq_trim_47[combat_seq_trim_47$Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.<0.05,]
combat_seq_p_47 <- combat_seq_p_47[is.finite(combat_seq_p_47$Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.),]
combat_seq_p_48 <- combat_seq_trim_48[combat_seq_trim_48$Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.<0.05,]
combat_seq_p_48 <- combat_seq_p_48[is.finite(combat_seq_p_48$Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.),]
combat_seq_p_137 <- combat_seq_trim_137[combat_seq_trim_137$Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.<0.05,]
combat_seq_p_137 <- combat_seq_p_137[is.finite(combat_seq_p_137$Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.),]
combat_seq_p_168 <- combat_seq_trim_168[combat_seq_trim_168$Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.<0.05,]
combat_seq_p_168 <- combat_seq_p_168[is.finite(combat_seq_p_168$Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.),]
combat_seq_p_173 <- combat_seq_trim_173[combat_seq_trim_173$Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.<0.05,]
combat_seq_p_173 <- combat_seq_p_173[is.finite(combat_seq_p_173$Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.),]
combat_seq_p_242 <- combat_seq_trim_242[combat_seq_trim_242$Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.<0.05,]
combat_seq_p_242 <- combat_seq_p_242[is.finite(combat_seq_p_242$Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.),]
combat_seq_p_245 <- combat_seq_trim_245[combat_seq_trim_245$Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.<0.05,]
combat_seq_p_245 <- combat_seq_p_245[is.finite(combat_seq_p_245$Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.),]

##
## Waiting to do DESeq2
##

# MBatch - EB
mbatch_eb_p_47 <- mbatch_eb_trim_47[mbatch_eb_trim_47$Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.<0.05,]
mbatch_eb_p_47 <- mbatch_eb_p_47[is.finite(mbatch_eb_p_47$Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.),]
mbatch_eb_p_48 <- mbatch_eb_trim_48[mbatch_eb_trim_48$Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.<0.05,]
mbatch_eb_p_48 <- mbatch_eb_p_48[is.finite(mbatch_eb_p_48$Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.),]
mbatch_eb_p_137 <- mbatch_eb_trim_137[mbatch_eb_trim_137$Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.<0.05,]
mbatch_eb_p_137 <- mbatch_eb_p_137[is.finite(mbatch_eb_p_137$Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.),]
mbatch_eb_p_168 <- mbatch_eb_trim_168[mbatch_eb_trim_168$Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.<0.05,]
mbatch_eb_p_168 <- mbatch_eb_p_168[is.finite(mbatch_eb_p_168$Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.),]
mbatch_eb_p_173 <- mbatch_eb_trim_173[mbatch_eb_trim_173$Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.<0.05,]
mbatch_eb_p_173 <- mbatch_eb_p_173[is.finite(mbatch_eb_p_173$Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.),]
mbatch_eb_p_242 <- mbatch_eb_trim_242[mbatch_eb_trim_242$Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.<0.05,]
mbatch_eb_p_242 <- mbatch_eb_p_242[is.finite(mbatch_eb_p_242$Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.),]
mbatch_eb_p_245 <- mbatch_eb_trim_245[mbatch_eb_trim_245$Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.<0.05,]
mbatch_eb_p_245 <- mbatch_eb_p_245[is.finite(mbatch_eb_p_245$Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.),]

# MBatch - AN
mbatch_an_p_47 <- mbatch_an_trim_47[mbatch_an_trim_47$Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.<0.05,]
mbatch_an_p_47 <- mbatch_an_p_47[is.finite(mbatch_an_p_47$Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.),]
mbatch_an_p_48 <- mbatch_an_trim_48[mbatch_an_trim_48$Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.<0.05,]
mbatch_an_p_48 <- mbatch_an_p_48[is.finite(mbatch_an_p_48$Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.),]
mbatch_an_p_137 <- mbatch_an_trim_137[mbatch_an_trim_137$Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.<0.05,]
mbatch_an_p_137 <- mbatch_an_p_137[is.finite(mbatch_an_p_137$Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.),]
mbatch_an_p_168 <- mbatch_an_trim_168[mbatch_an_trim_168$Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.<0.05,]
mbatch_an_p_168 <- mbatch_an_p_168[is.finite(mbatch_an_p_168$Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.),]
mbatch_an_p_173 <- mbatch_an_trim_173[mbatch_an_trim_173$Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.<0.05,]
mbatch_an_p_173 <- mbatch_an_p_173[is.finite(mbatch_an_p_173$Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.),]
mbatch_an_p_242 <- mbatch_an_trim_242[mbatch_an_trim_242$Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.<0.05,]
mbatch_an_p_242 <- mbatch_an_p_242[is.finite(mbatch_an_p_242$Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.),]
mbatch_an_p_245 <- mbatch_an_trim_245[mbatch_an_trim_245$Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.<0.05,]
mbatch_an_p_245 <- mbatch_an_p_245[is.finite(mbatch_an_p_245$Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.),]

# MBatch - MP
mbatch_mp_p_47 <- mbatch_mp_trim_47[mbatch_mp_trim_47$Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.<0.05,]
mbatch_mp_p_47 <- mbatch_mp_p_47[is.finite(mbatch_mp_p_47$Adj.p.value_.FLT...GLDS.47.v.GC...GLDS.47.),]
mbatch_mp_p_48 <- mbatch_mp_trim_48[mbatch_mp_trim_48$Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.<0.05,]
mbatch_mp_p_48 <- mbatch_mp_p_48[is.finite(mbatch_mp_p_48$Adj.p.value_.FLT...GLDS.48.v.GC...GLDS.48.),]
mbatch_mp_p_137 <- mbatch_mp_trim_137[mbatch_mp_trim_137$Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.<0.05,]
mbatch_mp_p_137 <- mbatch_mp_p_137[is.finite(mbatch_mp_p_137$Adj.p.value_.FLT...GLDS.137.v.GC...GLDS.137.),]
mbatch_mp_p_168 <- mbatch_mp_trim_168[mbatch_mp_trim_168$Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.<0.05,]
mbatch_mp_p_168 <- mbatch_mp_p_168[is.finite(mbatch_mp_p_168$Adj.p.value_.FLT...GLDS.168.v.GC...GLDS.168.),]
mbatch_mp_p_173 <- mbatch_mp_trim_173[mbatch_mp_trim_173$Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.<0.05,]
mbatch_mp_p_173 <- mbatch_mp_p_173[is.finite(mbatch_mp_p_173$Adj.p.value_.FLT...GLDS.173.v.GC...GLDS.173.),]
mbatch_mp_p_242 <- mbatch_mp_trim_242[mbatch_mp_trim_242$Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.<0.05,]
mbatch_mp_p_242 <- mbatch_mp_p_242[is.finite(mbatch_mp_p_242$Adj.p.value_.FLT...GLDS.242.v.GC...GLDS.242.),]
mbatch_mp_p_245 <- mbatch_mp_trim_245[mbatch_mp_trim_245$Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.<0.05,]
mbatch_mp_p_245 <- mbatch_mp_p_245[is.finite(mbatch_mp_p_245$Adj.p.value_.FLT...GLDS.245.v.GC...GLDS.245.),]


## Get total counts for all filtered DEGs from all datasets
## Determine how many match the final table
# prepare table headings/names
tot_row_names <- c("Original", "ComBat", "ComBat_seq", "MBatch-EB", "MBatch-AN", "MBatch-MP")
match_row_names <- c("ComBat", "ComBat_seq", "MBatch-EB", "MBatch-AN", "MBatch-MP")

# GLDS-47 comparison
orig_47_deg_count <- length(rownames(p_47)) # 23
combat_47_deg_count <- length(rownames(combat_p_47)) # 248
combat_seq_47_deg_count <- length(rownames(combat_seq_p_47))
# deseq2_47_deg_count <- length(rownames(deseq2_p_47))
mbatch_eb_47_deg_count <- length(rownames(mbatch_eb_p_47))
mbatch_an_47_deg_count <- length(rownames(mbatch_an_p_47))
mbatch_mp_47_deg_count <- length(rownames(mbatch_mp_p_47))
match_combat_47 <- sum(countMatches(rownames(p_47), rownames(combat_p_47))) # 9
match_combat_seq_47 <- sum(countMatches(rownames(p_47), rownames(combat_seq_p_47)))
# match_deseq2_47 <- sum(countMatches(rownames(p_47), rownames(deseq2_p_47)))
match_mbatch_eb_47 <- sum(countMatches(rownames(p_47), rownames(mbatch_eb_p_47)))
match_mbatch_an_47 <- sum(countMatches(rownames(p_47), rownames(mbatch_an_p_47)))
match_mbatch_mp_47 <- sum(countMatches(rownames(p_47), rownames(mbatch_mp_p_47)))

# remember to add deseq2 to these
tot_col_47 <- rbind(orig_47_deg_count, combat_47_deg_count, combat_seq_47_deg_count, mbatch_eb_47_deg_count, mbatch_an_47_deg_count, mbatch_mp_47_deg_count)
rownames(tot_col_47) <- tot_row_names
colnames(tot_col_47) <- "GLDS-47"

match_col_47 <- rbind(match_combat_47, match_combat_seq_47, match_mbatch_eb_47, match_mbatch_an_47, match_mbatch_mp_47)
rownames(match_col_47) <- match_row_names
colnames(match_col_47) <- "GLDS-47"


# GLDS-48 comparison
orig_48_deg_count <- length(rownames(p_48))
combat_48_deg_count <- length(rownames(combat_p_48))
combat_seq_48_deg_count <- length(rownames(combat_seq_p_48))
# deseq2_47_deg_count <- length(rownames(deseq2_p_47))
mbatch_eb_48_deg_count <- length(rownames(mbatch_eb_p_48))
mbatch_an_48_deg_count <- length(rownames(mbatch_an_p_48))
mbatch_mp_48_deg_count <- length(rownames(mbatch_mp_p_48))
match_combat_48 <- sum(countMatches(rownames(p_48), rownames(combat_p_48)))
match_combat_seq_48 <- sum(countMatches(rownames(p_48), rownames(combat_seq_p_48)))
# match_deseq2_47 <- sum(countMatches(rownames(p_47), rownames(deseq2_p_47)))
match_mbatch_eb_48 <- sum(countMatches(rownames(p_48), rownames(mbatch_eb_p_48)))
match_mbatch_an_48 <- sum(countMatches(rownames(p_48), rownames(mbatch_an_p_48)))
match_mbatch_mp_48 <- sum(countMatches(rownames(p_48), rownames(mbatch_mp_p_48)))

# remember to add deseq2 to these
tot_col_48 <- rbind(orig_48_deg_count, combat_48_deg_count, combat_seq_48_deg_count, mbatch_eb_48_deg_count, mbatch_an_48_deg_count, mbatch_mp_48_deg_count)
rownames(tot_col_48) <- tot_row_names
colnames(tot_col_48) <- "GLDS-48"

match_col_48 <- rbind(match_combat_48, match_combat_seq_48, match_mbatch_eb_48, match_mbatch_an_48, match_mbatch_mp_48)
rownames(match_col_48) <- match_row_names
colnames(match_col_48) <- "GLDS-48"


# # GLDS-48-C comparison
# orig_48_C_deg_count <- length(rownames(p_48_C))
# combat_48_C_deg_count <- length(rownames(combat_p_48))
# combat_seq_48_C_deg_count <- length(rownames(combat_seq_p_48))
# # deseq2_47_deg_count <- length(rownames(deseq2_p_47))
# mbatch_eb_48_C_deg_count <- length(rownames(mbatch_eb_p_48))
# mbatch_an_48_C_deg_count <- length(rownames(mbatch_an_p_48))
# mbatch_mp_48_C_deg_count <- length(rownames(mbatch_mp_p_48))
# match_combat_48_C <- sum(countMatches(rownames(p_48_C), rownames(combat_p_48)))
# match_combat_seq_48_C <- sum(countMatches(rownames(p_48_C), rownames(combat_seq_p_48)))
# # match_deseq2_47 <- sum(countMatches(rownames(p_47), rownames(deseq2_p_47)))
# match_mbatch_eb_48_C <- sum(countMatches(rownames(p_48_C), rownames(mbatch_eb_p_48)))
# match_mbatch_an_48_C <- sum(countMatches(rownames(p_48_C), rownames(mbatch_an_p_48)))
# match_mbatch_mp_48_C <- sum(countMatches(rownames(p_48_C), rownames(mbatch_mp_p_48)))
# 
# # remember to add deseq2 to these
# tot_col_48_C <- rbind(orig_48_C_deg_count, combat_48_C_deg_count, combat_seq_48_C_deg_count, mbatch_eb_48_C_deg_count, mbatch_an_48_C_deg_count, mbatch_mp_48_C_deg_count)
# rownames(tot_col_48_C) <- tot_row_names
# colnames(tot_col_48_C) <- "GLDS-48-C"
# 
# match_col_48_C <- rbind(match_combat_48_C, match_combat_seq_48_C, match_mbatch_eb_48_C, match_mbatch_an_48_C, match_mbatch_mp_48_C)
# rownames(match_col_48_C) <- match_row_names
# colnames(match_col_48_C) <- "GLDS-48-C"
# 
# 
# # GLDS-48-I comparison
# orig_48_I_deg_count <- length(rownames(p_48_I))
# combat_48_I_deg_count <- length(rownames(combat_p_48))
# combat_seq_48_I_deg_count <- length(rownames(combat_seq_p_48))
# # deseq2_47_deg_count <- length(rownames(deseq2_p_47))
# mbatch_eb_48_I_deg_count <- length(rownames(mbatch_eb_p_48))
# mbatch_an_48_I_deg_count <- length(rownames(mbatch_an_p_48))
# mbatch_mp_48_I_deg_count <- length(rownames(mbatch_mp_p_48))
# match_combat_48_I <- sum(countMatches(rownames(p_48_I), rownames(combat_p_48)))
# match_combat_seq_48_I <- sum(countMatches(rownames(p_48_I), rownames(combat_seq_p_48)))
# # match_deseq2_47 <- sum(countMatches(rownames(p_47), rownames(deseq2_p_47)))
# match_mbatch_eb_48_I <- sum(countMatches(rownames(p_48_I), rownames(mbatch_eb_p_48)))
# match_mbatch_an_48_I <- sum(countMatches(rownames(p_48_I), rownames(mbatch_an_p_48)))
# match_mbatch_mp_48_I <- sum(countMatches(rownames(p_48_I), rownames(mbatch_mp_p_48)))
# 
# # remember to add deseq2 to these
# tot_col_48_I <- rbind(orig_48_I_deg_count, combat_48_I_deg_count, combat_seq_48_I_deg_count, mbatch_eb_48_I_deg_count, mbatch_an_48_I_deg_count, mbatch_mp_48_I_deg_count)
# rownames(tot_col_48_I) <- tot_row_names
# colnames(tot_col_48_I) <- "GLDS-48-I"
# 
# match_col_48_I <- rbind(match_combat_48_I, match_combat_seq_48_I, match_mbatch_eb_48_I, match_mbatch_an_48_I, match_mbatch_mp_48_I)
# rownames(match_col_48_I) <- match_row_names
# colnames(match_col_48_I) <- "GLDS-48-I"


# GLDS-137 comparison
orig_137_deg_count <- length(rownames(p_137))
combat_137_deg_count <- length(rownames(combat_p_137))
combat_seq_137_deg_count <- length(rownames(combat_seq_p_137))
# deseq2_47_deg_count <- length(rownames(deseq2_p_47))
mbatch_eb_137_deg_count <- length(rownames(mbatch_eb_p_137))
mbatch_an_137_deg_count <- length(rownames(mbatch_an_p_137))
mbatch_mp_137_deg_count <- length(rownames(mbatch_mp_p_137))
match_combat_137 <- sum(countMatches(rownames(p_137), rownames(combat_p_137)))
match_combat_seq_137 <- sum(countMatches(rownames(p_137), rownames(combat_seq_p_137)))
# match_deseq2_47 <- sum(countMatches(rownames(p_47), rownames(deseq2_p_47)))
match_mbatch_eb_137 <- sum(countMatches(rownames(p_137), rownames(mbatch_eb_p_137)))
match_mbatch_an_137 <- sum(countMatches(rownames(p_137), rownames(mbatch_an_p_137)))
match_mbatch_mp_137 <- sum(countMatches(rownames(p_137), rownames(mbatch_mp_p_137)))

# remember to add deseq2 to these
tot_col_137 <- rbind(orig_137_deg_count, combat_137_deg_count, combat_seq_137_deg_count, mbatch_eb_137_deg_count, mbatch_an_137_deg_count, mbatch_mp_137_deg_count)
rownames(tot_col_137) <- tot_row_names
colnames(tot_col_137) <- "GLDS-137"

match_col_137 <- rbind(match_combat_137, match_combat_seq_137, match_mbatch_eb_137, match_mbatch_an_137, match_mbatch_mp_137)
rownames(match_col_137) <- match_row_names
colnames(match_col_137) <- "GLDS-137"


# GLDS-168 comparison
orig_168_deg_count <- length(rownames(p_168))
combat_168_deg_count <- length(rownames(combat_p_168))
combat_seq_168_deg_count <- length(rownames(combat_seq_p_168))
# deseq2_47_deg_count <- length(rownames(deseq2_p_47))
mbatch_eb_168_deg_count <- length(rownames(mbatch_eb_p_168))
mbatch_an_168_deg_count <- length(rownames(mbatch_an_p_168))
mbatch_mp_168_deg_count <- length(rownames(mbatch_mp_p_168))
match_combat_168 <- sum(countMatches(rownames(p_168), rownames(combat_p_168)))
match_combat_seq_168 <- sum(countMatches(rownames(p_168), rownames(combat_seq_p_168)))
# match_deseq2_47 <- sum(countMatches(rownames(p_47), rownames(deseq2_p_47)))
match_mbatch_eb_168 <- sum(countMatches(rownames(p_168), rownames(mbatch_eb_p_168)))
match_mbatch_an_168 <- sum(countMatches(rownames(p_168), rownames(mbatch_an_p_168)))
match_mbatch_mp_168 <- sum(countMatches(rownames(p_168), rownames(mbatch_mp_p_168)))

# remember to add deseq2 to these
tot_col_168 <- rbind(orig_168_deg_count, combat_168_deg_count, combat_seq_168_deg_count, mbatch_eb_168_deg_count, mbatch_an_168_deg_count, mbatch_mp_168_deg_count)
rownames(tot_col_168) <- tot_row_names
colnames(tot_col_168) <- "GLDS-168"

match_col_168 <- rbind(match_combat_168, match_combat_seq_168, match_mbatch_eb_168, match_mbatch_an_168, match_mbatch_mp_168)
rownames(match_col_168) <- match_row_names
colnames(match_col_168) <- "GLDS-168"


# # GLDS-168_1 comparison
# orig_168_1_deg_count <- length(rownames(p_168_1))
# combat_168_1_deg_count <- length(rownames(combat_p_168))
# combat_seq_168_1_deg_count <- length(rownames(combat_seq_p_168))
# # deseq2_47_deg_count <- length(rownames(deseq2_p_47))
# mbatch_eb_168_1_deg_count <- length(rownames(mbatch_eb_p_168))
# mbatch_an_168_1_deg_count <- length(rownames(mbatch_an_p_168))
# mbatch_mp_168_1_deg_count <- length(rownames(mbatch_mp_p_168))
# match_combat_168_1 <- sum(countMatches(rownames(p_168_1), rownames(combat_p_168)))
# match_combat_seq_168_1 <- sum(countMatches(rownames(p_168_1), rownames(combat_seq_p_168)))
# # match_deseq2_47 <- sum(countMatches(rownames(p_47), rownames(deseq2_p_47)))
# match_mbatch_eb_168_1 <- sum(countMatches(rownames(p_168_1), rownames(mbatch_eb_p_168)))
# match_mbatch_an_168_1 <- sum(countMatches(rownames(p_168_1), rownames(mbatch_an_p_168)))
# match_mbatch_mp_168_1 <- sum(countMatches(rownames(p_168_1), rownames(mbatch_mp_p_168)))
# 
# # remember to add deseq2 to these
# tot_col_168_1 <- rbind(orig_168_1_deg_count, combat_168_1_deg_count, combat_seq_168_1_deg_count, mbatch_eb_168_1_deg_count, mbatch_an_168_1_deg_count, mbatch_mp_168_1_deg_count)
# rownames(tot_col_168_1) <- tot_row_names
# colnames(tot_col_168_1) <- "GLDS-168-noERCC"
# 
# match_col_168_1 <- rbind(match_combat_168_1, match_combat_seq_168_1, match_mbatch_eb_168_1, match_mbatch_an_168_1, match_mbatch_mp_168_1)
# rownames(match_col_168_1) <- match_row_names
# colnames(match_col_168_1) <- "GLDS-168-noERCC"
# 
# 
# # GLDS-168_3 comparison
# orig_168_3_deg_count <- length(rownames(p_168_3))
# combat_168_3_deg_count <- length(rownames(combat_p_168))
# combat_seq_168_3_deg_count <- length(rownames(combat_seq_p_168))
# # deseq2_47_deg_count <- length(rownames(deseq2_p_47))
# mbatch_eb_168_3_deg_count <- length(rownames(mbatch_eb_p_168))
# mbatch_an_168_3_deg_count <- length(rownames(mbatch_an_p_168))
# mbatch_mp_168_3_deg_count <- length(rownames(mbatch_mp_p_168))
# match_combat_168_3 <- sum(countMatches(rownames(p_168_3), rownames(combat_p_168)))
# match_combat_seq_168_3 <- sum(countMatches(rownames(p_168_3), rownames(combat_seq_p_168)))
# # match_deseq2_47 <- sum(countMatches(rownames(p_47), rownames(deseq2_p_47)))
# match_mbatch_eb_168_3 <- sum(countMatches(rownames(p_168_3), rownames(mbatch_eb_p_168)))
# match_mbatch_an_168_3 <- sum(countMatches(rownames(p_168_3), rownames(mbatch_an_p_168)))
# match_mbatch_mp_168_3 <- sum(countMatches(rownames(p_168_3), rownames(mbatch_mp_p_168)))
# 
# # remember to add deseq2 to these
# tot_col_168_3 <- rbind(orig_168_3_deg_count, combat_168_3_deg_count, combat_seq_168_3_deg_count, mbatch_eb_168_3_deg_count, mbatch_an_168_3_deg_count, mbatch_mp_168_3_deg_count)
# rownames(tot_col_168_3) <- tot_row_names
# colnames(tot_col_168_3) <- "GLDS-168-wERCC"
# 
# match_col_168_3 <- rbind(match_combat_168_3, match_combat_seq_168_3, match_mbatch_eb_168_3, match_mbatch_an_168_3, match_mbatch_mp_168_3)
# rownames(match_col_168_3) <- match_row_names
# colnames(match_col_168_3) <- "GLDS-168-wERCC"


# GLDS-173 comparison
orig_173_deg_count <- length(rownames(p_173))
combat_173_deg_count <- length(rownames(combat_p_173))
combat_seq_173_deg_count <- length(rownames(combat_seq_p_173))
# deseq2_47_deg_count <- length(rownames(deseq2_p_47))
mbatch_eb_173_deg_count <- length(rownames(mbatch_eb_p_173))
mbatch_an_173_deg_count <- length(rownames(mbatch_an_p_173))
mbatch_mp_173_deg_count <- length(rownames(mbatch_mp_p_173))
match_combat_173 <- sum(countMatches(rownames(p_173), rownames(combat_p_173)))
match_combat_seq_173 <- sum(countMatches(rownames(p_173), rownames(combat_seq_p_173)))
# match_deseq2_47 <- sum(countMatches(rownames(p_47), rownames(deseq2_p_47)))
match_mbatch_eb_173 <- sum(countMatches(rownames(p_173), rownames(mbatch_eb_p_173)))
match_mbatch_an_173 <- sum(countMatches(rownames(p_173), rownames(mbatch_an_p_173)))
match_mbatch_mp_173 <- sum(countMatches(rownames(p_173), rownames(mbatch_mp_p_173)))

# remember to add deseq2 to these
tot_col_173 <- rbind(orig_173_deg_count, combat_173_deg_count, combat_seq_173_deg_count, mbatch_eb_173_deg_count, mbatch_an_173_deg_count, mbatch_mp_173_deg_count)
rownames(tot_col_173) <- tot_row_names
colnames(tot_col_173) <- "GLDS-173"

match_col_173 <- rbind(match_combat_173, match_combat_seq_173, match_mbatch_eb_173, match_mbatch_an_173, match_mbatch_mp_173)
rownames(match_col_173) <- match_row_names
colnames(match_col_173) <- "GLDS-173"


# GLDS-242 comparison
orig_242_deg_count <- length(rownames(p_242))
combat_242_deg_count <- length(rownames(combat_p_242))
combat_seq_242_deg_count <- length(rownames(combat_seq_p_242))
# deseq2_47_deg_count <- length(rownames(deseq2_p_47))
mbatch_eb_242_deg_count <- length(rownames(mbatch_eb_p_242))
mbatch_an_242_deg_count <- length(rownames(mbatch_an_p_242))
mbatch_mp_242_deg_count <- length(rownames(mbatch_mp_p_242))
match_combat_242 <- sum(countMatches(rownames(p_242), rownames(combat_p_242)))
match_combat_seq_242 <- sum(countMatches(rownames(p_242), rownames(combat_seq_p_242)))
# match_deseq2_47 <- sum(countMatches(rownames(p_47), rownames(deseq2_p_47)))
match_mbatch_eb_242 <- sum(countMatches(rownames(p_242), rownames(mbatch_eb_p_242)))
match_mbatch_an_242 <- sum(countMatches(rownames(p_242), rownames(mbatch_an_p_242)))
match_mbatch_mp_242 <- sum(countMatches(rownames(p_242), rownames(mbatch_mp_p_242)))

# remember to add deseq2 to these
tot_col_242 <- rbind(orig_242_deg_count, combat_242_deg_count, combat_seq_242_deg_count, mbatch_eb_242_deg_count, mbatch_an_242_deg_count, mbatch_mp_242_deg_count)
rownames(tot_col_242) <- tot_row_names
colnames(tot_col_242) <- "GLDS-242"

match_col_242 <- rbind(match_combat_242, match_combat_seq_242, match_mbatch_eb_242, match_mbatch_an_242, match_mbatch_mp_242)
rownames(match_col_242) <- match_row_names
colnames(match_col_242) <- "GLDS-242"



# GLDS-245 comparison
orig_245_deg_count <- length(rownames(p_245))
combat_245_deg_count <- length(rownames(combat_p_245))
combat_seq_245_deg_count <- length(rownames(combat_seq_p_245))
# deseq2_47_deg_count <- length(rownames(deseq2_p_47))
mbatch_eb_245_deg_count <- length(rownames(mbatch_eb_p_245))
mbatch_an_245_deg_count <- length(rownames(mbatch_an_p_245))
mbatch_mp_245_deg_count <- length(rownames(mbatch_mp_p_245))
match_combat_245 <- sum(countMatches(rownames(p_245), rownames(combat_p_245)))
match_combat_seq_245 <- sum(countMatches(rownames(p_245), rownames(combat_seq_p_245)))
# match_deseq2_47 <- sum(countMatches(rownames(p_47), rownames(deseq2_p_47)))
match_mbatch_eb_245 <- sum(countMatches(rownames(p_245), rownames(mbatch_eb_p_245)))
match_mbatch_an_245 <- sum(countMatches(rownames(p_245), rownames(mbatch_an_p_245)))
match_mbatch_mp_245 <- sum(countMatches(rownames(p_245), rownames(mbatch_mp_p_245)))

# remember to add deseq2 to these
tot_col_245 <- rbind(orig_245_deg_count, combat_245_deg_count, combat_seq_245_deg_count, mbatch_eb_245_deg_count, mbatch_an_245_deg_count, mbatch_mp_245_deg_count)
rownames(tot_col_245) <- tot_row_names
colnames(tot_col_245) <- "GLDS-245"

match_col_245 <- rbind(match_combat_245, match_combat_seq_245, match_mbatch_eb_245, match_mbatch_an_245, match_mbatch_mp_245)
rownames(match_col_245) <- match_row_names
colnames(match_col_245) <- "GLDS-245"

## Combine tables
total_table <- cbind(tot_col_47, tot_col_48, tot_col_137, tot_col_168, tot_col_173, tot_col_242, tot_col_245)
match_table <- cbind(match_col_47, match_col_48, match_col_137, match_col_168, match_col_173, match_col_242, match_col_245)

total_no_orig_table <- total_table[2:6,]
extra_table <- total_no_orig_table - match_table


## total_table contains count of total DEGs
## match_table contains count of DEGs that match those of original dataset
## extra_table contains count of DEGs that don't match those of original


