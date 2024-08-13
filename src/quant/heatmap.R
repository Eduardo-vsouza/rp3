# library("DESeq2")
# library("RColorBrewer")
# library("gplots")
# library("pheatmap")
# library("dplyr")
# library('glue')
# library('data.table')
# library(dplyr)
#
# # library(reshape2)
#
# args <- commandArgs(trailingOnly = TRUE)
# # args: <counts> <outdir> <correlation_cutoff> <padj_cutoff>
# #cor_table <- as.data.frame(read.table(args[1], sep='\t'))
# counts_file <- args[1]
# samples_file <- args[2]
# outdir <- args[3]
# # library("DESeq2")
#
# # library('data.table')
# # library(dplyr)
# library(ComplexHeatmap)
# library(circlize)
#
#
# diff <- as.matrix(read.csv(counts_file, sep='\t', row.names=1))
# coldata <- read.csv(samples_file, sep='\t', row.names=1)
# annocol <- coldata[!coldata$Biorep]
# annocol <- coldata[!coldata$Fraction]
# annocol <- coldata[!coldata$Techrep]
# annocol
# ha <- HeatmapAnnotation(
#   Condition = coldata$Condition,
#   col = list(Condition = c("UI" = "red", '1640' = 'green')))
#
# pdf(paste0(outdir, "_heatmap.pdf"), width = 8, height = 6)
#
# Heatmap(
#   diff,
#   name = "Expression",
#   top_annotation = ha,
#   show_column_names = TRUE,
#   show_row_names = TRUE
# )
# dev.off()
#


##############################



# library("DESeq2")
# library("RColorBrewer")
# library("gplots")
# library("pheatmap")
# library("dplyr")
# library('glue')
# library('data.table')
# library(dplyr)
#
# library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
# args: <counts> <outdir> <correlation_cutoff> <padj_cutoff>
#cor_table <- as.data.frame(read.table(args[1], sep='\t'))
counts_file <- args[1]
samples_file <- args[2]
outdir <- args[3]
group1 <- args[4]
group2 <- args[5]
# library("DESeq2")
# counts_file <-  "/home/microway/projects/serva_IP_e4orf6-7/rp3_IP_comparisons_genomeWide/240805_newSamples_oldSamples_flashLFQ/flashLFQ/comparisons/1640_x_102/standard_intensities_heatmap.csv"
# samples_file <- '/home/microway/projects/serva_IP_e4orf6-7/rp3_IP_comparisons_genomeWide/240805_newSamples_oldSamples_flashLFQ/flashLFQ/cat_mzml_files/ExperimentalDesign.tsv'
# library('data.table')
# library(dplyr)
library(ComplexHeatmap)
library(circlize)
colnames(diff)

diff <- as.matrix(read.csv(counts_file, sep='\t', row.names=1))
coldata <- read.csv(samples_file, sep='\t', row.names=1)
coldata
# group1 = group1
# group2 = group2
colnames(diff)
colnames(diff) <- gsub("^X", "", colnames(diff))

coldata <- subset(coldata, Condition %in% c(group1, group2))
coldata
ha <- HeatmapAnnotation(
  Condition = coldata$Condition,
  col = list(Condition = setNames(c("red", "blue"), c(group1, group2)))
)

pdf(paste0(outdir, "_heatmap.pdf"), width = 8, height = 6)

Heatmap(
  diff,
  name = "Expression",
  top_annotation = ha,
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),

)

dev.off()

