#tell R where the libraries are
.libPaths("/mnt/beegfs/marcovicari3/Rlib")
library(devtools)
library(ggplot2)
library(cowplot)
library(Matrix)
library(MatrixModels)
library(plyr)
library(dplyr)
library(Seurat)
library(monocle)
library(reshape2)
library(stringr)
library(netbiov)

rm(list = ls()) # clear the environment 
options(warn=-1) # turn off warning message globally 
options(stringsAsFactors = FALSE)

# set working directory ####
setwd("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/")

# load objects ####
load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/AL1_seurat_filt_scale")
load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/AL2_seurat_filt_scale")
load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/AL3_seurat_filt_scale")
load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/CPN2_seurat_filt_scale")
load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/CPN3_seurat_filt_scale")
load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/CPN4_seurat_filt_scale")
load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/CPP1_seurat_filt_scale")
load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/CPP3_seurat_filt_scale")
load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/CPP4_seurat_filt_scale")
load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/CPP.combined_seurat_filt_merge_scale_batch")
load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/CPN.combined_seurat_filt_merge_scale_batch")
load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/AL.combined_seurat_filt_merge_scale_batch")
# load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/HAM.combined_unaligned")
load("/mnt/beegfs/marcovicari3/singlecell/10xsc/seurat/objects/HAM_agg_merge_cca_aligned")

# ### Acquire data ####
#  
# CPN2_ <-("/mnt/beegfs/marcovicari3/singlecell/10xsc/data/170629_NB501241_0017_AH253VBGX3_CPN2/2_count_outs/outs/filtered_gene_bc_matrices/GRCh38/")
# Data2 <-Read10X(data.dir=CPN2_)
# colnames(Data2) <- paste("CPN2", colnames(Data2), sep = "_")
# CPN2 <- CreateSeuratObject(raw.data = Data2, min.cells = 5)
# # CPN2@data <- CPN2@data[,sample(1:ncol(CPN2@data),200)]
# CPN2@meta.data$sample<-"CPN2"
# CPN2@meta.data$batch<-"2"
# CPN2@meta.data$condition<-"CPN"
# 
# MyDirCPN3<-"/mnt/beegfs/marcovicari3/singlecell/10xsc/data/CPN3/2_count_outs/outs/filtered_gene_bc_matrices/GRCh38/"
# Data9<-Read10X(data.dir=MyDirCPN3)
# colnames(Data9) <- paste("CPN3", colnames(Data9), sep = "_")
# CPN3 <- CreateSeuratObject(raw.data = Data9, min.cells = 5)
# # CPN3@data <- CPN3@data[,sample(1:ncol(CPN3@data),200)]
# CPN3@meta.data$sample<-"CPN3"
# CPN3@meta.data$batch<-"3"
# CPN3@meta.data$condition<-"CPN"
# 
# MyDirCPN4<-("/mnt/beegfs/marcovicari3/singlecell/10xsc/data/CPN4/2_count_outs/outs/filtered_gene_bc_matrices/GRCh38/")
# Data10<-Read10X(data.dir=MyDirCPN4)
# colnames(Data10) <- paste("CPN4", colnames(Data10), sep = "_")
# CPN4 <- CreateSeuratObject(raw.data = Data10, min.cells = 5)
# # CPN4@data <- CPN4@data[,sample(1:ncol(CPN4@data),200)]
# CPN4@meta.data$sample<-"CPN4"
# CPN4@meta.data$batch<-"1"
# CPN4@meta.data$condition<-"CPN"
# 
# MyDirAL1<-("/mnt/beegfs/marcovicari3/singlecell/10xsc/data/170718_NB501241_0024_AH37H2BGX3_AL/2_count_outs/outs/filtered_gene_bc_matrices/GRCh38/")
# Data3<-Read10X(data.dir=MyDirAL1)
# colnames(Data3) <- paste("AL1", colnames(Data3), sep = "_")
# AL1 <- CreateSeuratObject(raw.data = Data3, min.cells = 5)
# # AL1@data <- AL1@data[,sample(1:ncol(AL1@data),200)]
# AL1@meta.data$sample<-"AL1"
# AL1@meta.data$batch<-"1"
# AL1@meta.data$condition<-"AL"
# 
# MyDirAL2<-("/mnt/beegfs/marcovicari3/singlecell/10xsc/data/171211_NB501241_0065_AHYVHGBGX3_AL2/2_count_outs_B12/outs/filtered_gene_bc_matrices/GRCh38/")
# Data4<-Read10X(data.dir=MyDirAL2)
# colnames(Data4) <- paste("AL2", colnames(Data4), sep = "_")
# AL2 <- CreateSeuratObject(raw.data = Data4, min.cells = 5)
# # AL2@data <- AL2@data[,sample(1:ncol(AL2@data),200)]
# AL2@meta.data$sample<-"AL2"
# AL2@meta.data$batch<-"2"
# AL2@meta.data$condition<-"AL"
# 
# MyDirAL3<-("/mnt/beegfs/marcovicari3/singlecell/10xsc/data/171215_NB501241_0068_AHCHWCBGX5_AL3/2_count_outs/outs/filtered_gene_bc_matrices/GRCh38/")
# Data5<-Read10X(data.dir=MyDirAL3)
# colnames(Data5) <- paste("AL3", colnames(Data5), sep = "_")
# AL3 <- CreateSeuratObject(raw.data = Data5, min.cells = 5)
# # AL3@data <- AL3@data[,sample(1:ncol(AL3@data),200)]
# AL3@meta.data$sample<-"AL3"
# AL3@meta.data$batch<-"3"
# AL3@meta.data$condition<-"AL"
# 
# MyDirCPP1<-("/mnt/beegfs/marcovicari3/singlecell/10xsc/data/170626_NB501241_0016_AHYWF7BGX2_CPP/2_count_outs/outs/filtered_gene_bc_matrices/GRCh38/")
# Data6<-Read10X(data.dir=MyDirCPP1)
# colnames(Data6) <- paste("CPP1", colnames(Data6), sep = "_")
# CPP1 <- CreateSeuratObject(raw.data = Data6, min.cells = 5)
# # CPP1@data <- CPP1@data[,sample(1:ncol(CPP1@data),200)]
# CPP1@meta.data$sample<-"CPP1"
# CPP1@meta.data$batch<-"1"
# CPP1@meta.data$condition<-"CPP"
# 
# MyDirCPP3<-("/mnt/beegfs/marcovicari3/singlecell/10xsc/data/171214_NB501241_0067_AHYV2YBGX3_CPP3/2_count_outs/outs/filtered_gene_bc_matrices/GRCh38/")
# Data8<-Read10X(data.dir=MyDirCPP3)
# colnames(Data8) <- paste("CPP3", colnames(Data8), sep = "_")
# CPP3 <- CreateSeuratObject(raw.data = Data8, min.cells = 5)
# # CPP3@data <- CPP3@data[,sample(1:ncol(CPP3@data),200)]
# CPP3@meta.data$sample<-"CPP3"
# CPP3@meta.data$batch<-"3"
# CPP3@meta.data$condition<-"CPP"
# 
# MyDirCPP4<-("/mnt/beegfs/marcovicari3/singlecell/10xsc/data/CPP4/2_count_outs/outs/filtered_gene_bc_matrices/GRCh38/")
# Data11<-Read10X(data.dir=MyDirCPP4)
# colnames(Data11) <- paste("CPP4", colnames(Data11), sep = "_")
# CPP4 <- CreateSeuratObject(raw.data = Data11, min.cells = 5)
# # CPP4@data <- CPP4@data[,sample(1:ncol(CPP4@data),200)]
# CPP4@meta.data$sample<-"CPP4"
# CPP4@meta.data$batch<-"2"
# CPP4@meta.data$condition<-"CPP"
# 
# ### QC1. Calculate % mitochondrial genes ####
# # % of genes per cell, UMI per cell, mithocondrial genes and filtering of cells
# 
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = CPN2@data), value = TRUE)
# percent.mito <- Matrix::colSums(CPN2@raw.data[mito.genes, ])/Matrix::colSums(CPN2@raw.data)
# CPN2 <- AddMetaData(object = CPN2, metadata = percent.mito, col.name = "percent.mito")
# pdf("qc/vln/vln_CPN2_seurat.pdf")
# VlnPlot(object = CPN2, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
# dev.off()
# 
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = CPN3@data), value = TRUE)
# percent.mito <- Matrix::colSums(CPN3@raw.data[mito.genes, ])/Matrix::colSums(CPN3@raw.data)
# CPN3 <- AddMetaData(object = CPN3, metadata = percent.mito, col.name = "percent.mito")
# pdf("qc/vln/vln_CPN3_seurat.pdf")
# VlnPlot(object = CPN3, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
# dev.off()
# 
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = CPN4@data), value = TRUE)
# percent.mito <- Matrix::colSums(CPN4@raw.data[mito.genes, ])/Matrix::colSums(CPN4@raw.data)
# CPN4 <- AddMetaData(object = CPN4, metadata = percent.mito, col.name = "percent.mito")
# pdf("qc/vln/vln_CPN4_seurat.pdf")
# VlnPlot(object = CPN4, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
# dev.off()
# 
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = AL1@data), value = TRUE)
# percent.mito <- Matrix::colSums(AL1@raw.data[mito.genes, ])/Matrix::colSums(AL1@raw.data)
# AL1 <- AddMetaData(object = AL1, metadata = percent.mito, col.name = "percent.mito")
# pdf("qc/vln/vln_AL1_seurat.pdf")
# VlnPlot(object = AL1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
# dev.off()
# 
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = AL2@data), value = TRUE)
# percent.mito <- Matrix::colSums(AL2@raw.data[mito.genes, ])/Matrix::colSums(AL2@raw.data)
# AL2 <- AddMetaData(object = AL2, metadata = percent.mito, col.name = "percent.mito")
# pdf("qc/vln/vln_AL2_seurat.pdf")
# VlnPlot(object = AL2, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
# dev.off()
# 
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = AL3@data), value = TRUE)
# percent.mito <- Matrix::colSums(AL3@raw.data[mito.genes, ])/Matrix::colSums(AL3@raw.data)
# AL3 <- AddMetaData(object = AL3, metadata = percent.mito, col.name = "percent.mito")
# pdf("qc/vln/vln_AL3_agg_filt.pdf")
# VlnPlot(object = AL3, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
# dev.off()
# 
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = CPP1@data), value = TRUE)
# percent.mito <- Matrix::colSums(CPP1@raw.data[mito.genes, ])/Matrix::colSums(CPP1@raw.data)
# CPP1 <- AddMetaData(object = CPP1, metadata = percent.mito, col.name = "percent.mito")
# pdf("qc/vln/vln_CPP1_agg_filt.pdf")
# VlnPlot(object = CPP1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
# dev.off()
# 
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = CPP3@data), value = TRUE)
# percent.mito <- Matrix::colSums(CPP3@raw.data[mito.genes, ])/Matrix::colSums(CPP3@raw.data)
# CPP3 <- AddMetaData(object = CPP3, metadata = percent.mito, col.name = "percent.mito")
# pdf("qc/vln/vln_CPP3_agg_filt.pdf")
# VlnPlot(object = CPP3, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
# dev.off()
# 
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = CPP4@data), value = TRUE)
# percent.mito <- Matrix::colSums(CPP4@raw.data[mito.genes, ])/Matrix::colSums(CPP4@raw.data)
# CPP4 <- AddMetaData(object = CPP4, metadata = percent.mito, col.name = "percent.mito")
# pdf("qc/vln/vln_CPP4_agg_filt.pdf")
# VlnPlot(object = CPP4, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
# dev.off()
# 
# # QC2. Plot genes vs UMI and vs precent.mito ####
# 
# pdf("qc/geneplots/geneplot_CPN2_seurat.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = CPN2, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = CPN2, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_CPN3_seurat.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = CPN3, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = CPN3, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_CPN4_seurat.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = CPN4, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = CPN4, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_AL1_seurat.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = AL1, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = AL1, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_AL2_seurat.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = AL2, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = AL2, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_AL3_seurat.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = AL3, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = AL3, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_CPP1_seurat.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = CPP1, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = CPP1, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_CPP3_seurat.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = CPP3, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = CPP3, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_CPP4_seurat.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = CPP4, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = CPP4, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# ## save raw data ####
# saveRDS(CPN2, file="objects/CPN2_seurat.rds")
# saveRDS(CPN3, file="objects/CPN3_seurat.rds")
# saveRDS(CPN4, file="objects/CPN4_seurat.rds")
# saveRDS(CPP1, file="objects/CPP1_seurat.rds")
# saveRDS(CPP3, file="objects/CPP3_seurat.rds")
# saveRDS(CPP4, file="objects/CPP4_seurat.rds")
# saveRDS(AL1, file="objects/AL1_seurat.rds")
# saveRDS(AL2, file="objects/AL2_seurat.rds")
# saveRDS(AL3, file="objects/AL3_seurat.rds")
# 
# # Filter cells ####
# filteredAL1 <- FilterCells(object = AL1, subset.names = c("nGene","nUMI", "percent.mito"), low.thresholds = c(-250, -Inf, -Inf), high.thresholds = c(1500, 20000, 0.05))
# filteredAL2 <- FilterCells(object = AL2, subset.names = c("nGene","nUMI", "percent.mito"), low.thresholds = c(-250, -Inf, -Inf), high.thresholds = c(3250, 50000, 0.05))
# filteredAL3 <- FilterCells(object = AL3, subset.names = c("nGene","nUMI", "percent.mito"), low.thresholds = c(-250, -Inf, -Inf), high.thresholds = c(1750, 20000, 0.05))
# filteredCPN2 <- FilterCells(object = CPN2, subset.names = c("nGene","nUMI", "percent.mito"), low.thresholds = c(-300, -Inf, -Inf), high.thresholds = c(2000, 25000, 0.05))
# filteredCPN3 <- FilterCells(object = CPN3, subset.names = c("nGene","nUMI", "percent.mito"), low.thresholds = c(-250, -Inf, -Inf), high.thresholds = c(2000, 20000, 0.05))
# filteredCPN4 <- FilterCells(object = CPN4, subset.names = c("nGene","nUMI", "percent.mito"), low.thresholds = c(-300, -Inf, -Inf), high.thresholds = c(2000, 25000, 0.05))
# filteredCPP1 <- FilterCells(object = CPP1, subset.names = c("nGene","nUMI", "percent.mito"), low.thresholds = c(-250, -Inf, -Inf), high.thresholds = c(3000, 40000, 0.05))
# filteredCPP3 <- FilterCells(object = CPP3, subset.names = c("nGene","nUMI", "percent.mito"), low.thresholds = c(-100, -Inf, -Inf), high.thresholds = c(700, 3250, 0.05))
# filteredCPP4 <- FilterCells(object = CPP4, subset.names = c("nGene","nUMI", "percent.mito"), low.thresholds = c(-300, -Inf, -Inf), high.thresholds = c(2000, 12500, 0.05))
# 
# # QC3. Plot genes vs UMI and vs precent.mito post filtering ####
# 
# pdf("qc/geneplots/geneplot_CPN2_seurat_filt.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = CPN2, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = CPN2, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_CPN3_seurat_filt.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = CPN3, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = CPN3, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_CPN4_seurat_filt.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = CPN4, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = CPN4, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_AL1_seurat_filt.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = AL1, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = AL1, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_AL2_seurat_filt.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = AL2, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = AL2, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_AL3_seurat_filt.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = AL3, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = AL3, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_CPP1_seurat_filt.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = CPP1, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = CPP1, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_CPP3_seurat_filt.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = CPP3, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = CPP3, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# pdf("qc/geneplots/geneplot_CPP4_seurat_filt.pdf")
# par(mfrow = c(1, 2))
# GenePlot(object = CPP4, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = CPP4, gene1 = "nUMI", gene2 = "nGene")
# dev.off()
# 
# saveRDS(filteredCPN2, file="objects/CPN2_seurat_filt.rds")
# saveRDS(filteredCPN3, file="objects/CPN3_seurat_filt.rds")
# saveRDS(filteredCPN4, file="objects/CPN4_seurat_filt.rds")
# saveRDS(filteredCPP1, file="objects/CPP1_seurat_filt.rds")
# saveRDS(filteredCPP3, file="objects/CPP3_seurat_filt.rds")
# saveRDS(filteredCPP4, file="objects/CPP4_seurat_filt.rds")
# saveRDS(filteredAL1, file="objects/AL1_seurat_filt.rds")
# saveRDS(filteredAL2, file="objects/AL2_seurat_filt.rds")
# saveRDS(filteredAL3, file="objects/AL3_seurat_filt.rds")
# 
# 
# #Merge conditions
# AL.combined <- MergeSeurat(filteredAL1, filteredAL2, do.normalize = F)
# AL.combined <- MergeSeurat(AL.combined, filteredAL3, do.normalize = F)
# 
# CPN.combined <- MergeSeurat(filteredCPN2, filteredCPN3, do.normalize = F)
# CPN.combined <- MergeSeurat(CPN.combined, filteredCPN4, do.normalize = F)
# 
# CPP.combined <- MergeSeurat(filteredCPP1, filteredCPP3, do.normalize = F)
# CPP.combined <- MergeSeurat(CPP.combined, filteredCPP4, do.normalize = F)
# 
# # Normalize
# filteredCPN2 <- NormalizeData(object = filteredCPN2, normalization.method = "LogNormalize", scale.factor = 10000)
# filteredCPN3 <- NormalizeData(object = filteredCPN3, normalization.method = "LogNormalize", scale.factor = 10000)
# filteredCPN4 <- NormalizeData(object = filteredCPN4, normalization.method = "LogNormalize", scale.factor = 10000)
# filteredAL1 <- NormalizeData(object = filteredAL1, normalization.method = "LogNormalize", scale.factor = 10000)
# filteredAL2 <- NormalizeData(object = filteredAL2, normalization.method = "LogNormalize", scale.factor = 10000)
# filteredAL3 <- NormalizeData(object = filteredAL3, normalization.method = "LogNormalize", scale.factor = 10000)
# filteredCPP1 <- NormalizeData(object = filteredCPP1, normalization.method = "LogNormalize", scale.factor = 10000)
# filteredCPP3 <- NormalizeData(object = filteredCPP3, normalization.method = "LogNormalize", scale.factor = 10000)
# filteredCPP4 <- NormalizeData(object = filteredCPP4, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# AL.combined <- NormalizeData(object = AL.combined, normalization.method = "LogNormalize", scale.factor = 10000)
# CPN.combined <- NormalizeData(object = CPN.combined, normalization.method = "LogNormalize", scale.factor = 10000)
# CPP.combined <- NormalizeData(object = CPP.combined, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# # # Find variable genes
# # pdf("qc/variablegenes/variablegenes_CPN2_agg.pdf")
# filteredCPN2 <- FindVariableGenes(object = filteredCPN2, do.contour=FALSE)
# # dev.off()
# # pdf("qc/variablegenes/variablegenes_CPN3_agg.pdf")
# filteredCPN3 <- FindVariableGenes(object = filteredCPN3, do.contour=FALSE)
# # dev.off()
# # pdf("qc/variablegenes/variablegenes_CPN4_agg.pdf")
# filteredCPN4 <- FindVariableGenes(object = filteredCPN4, do.contour=FALSE)
# # dev.off()
# # pdf("qc/variablegenes/variablegenes_AL1_agg.pdf")
# filteredAL1 <- FindVariableGenes(object = filteredAL1, do.contour=FALSE)
# # dev.off()
# # pdf("qc/variablegenes/variablegenes_AL2_agg.pdf")
# filteredAL2 <- FindVariableGenes(object = filteredAL2, do.contour=FALSE)
# # dev.off()
# # pdf("qc/variablegenes/variablegenes_AL3_agg.pdf")
# filteredAL3 <- FindVariableGenes(object = filteredAL3, do.contour=FALSE)
# # dev.off()
# # pdf("qc/variablegenes/variablegenes_CPP1_agg.pdf")
# filteredCPP1 <- FindVariableGenes(object = filteredCPP1, do.contour=FALSE)
# # dev.off()
# # pdf("qc/variablegenes/variablegenes_CPP3_agg.pdf")
# filteredCPP3 <- FindVariableGenes(object = filteredCPP3, do.contour=FALSE)
# # dev.off()
# # pdf("qc/variablegenes/variablegenes_CPP4_agg.pdf")
# filteredCPP4 <- FindVariableGenes(object = filteredCPP4, do.contour=FALSE)
# # dev.off()
# # pdf("qc/variablegenes/variablegenes_AL_agg.pdf")
# AL.combined <- FindVariableGenes(object = AL.combined, do.contour=FALSE)
# # dev.off()
# # pdf("qc/variablegenes/variablegenes_CPN_agg.pdf")
# CPN.combined <- FindVariableGenes(object = CPN.combined, do.contour=FALSE)
# # dev.off()
# # pdf("qc/variablegenes/variablegenes_CPP_agg.pdf")
# CPP.combined <- FindVariableGenes(object = CPP.combined, do.contour=FALSE)
# # dev.off()

hvg.filteredCPN2<-rownames(head(filteredCPN2@hvg.info,2000))
hvg.filteredCPN3<-rownames(head(filteredCPN3@hvg.info,2000))
hvg.filteredCPN4<-rownames(head(filteredCPN4@hvg.info,2000))
hvg.filteredAL1<-rownames(head(filteredAL1@hvg.info,2000))
hvg.filteredAL2<-rownames(head(filteredAL2@hvg.info,2000))
hvg.filteredAL3<-rownames(head(filteredAL3@hvg.info,2000))
hvg.filteredCPP1<-rownames(head(filteredCPP1@hvg.info,2000))
hvg.filteredCPP3<-rownames(head(filteredCPP3@hvg.info,2000))
hvg.filteredCPP4<-rownames(head(filteredCPP4@hvg.info,2000))
hvg.AL.combined <-rownames(head(AL.combined@hvg.info,2000))
hvg.CPN.combined <-rownames(head(CPN.combined@hvg.info,2000))
hvg.CPP.combined <-rownames(head(CPP.combined@hvg.info,2000))

# filteredCPN2 <- ScaleData(object = filteredCPN2, vars.to.regress = c("nUMI", "percent.mito"))
# filteredCPN3 <- ScaleData(object = filteredCPN3, vars.to.regress = c("nUMI", "percent.mito"))
# filteredCPN4 <- ScaleData(object = filteredCPN4, vars.to.regress = c("nUMI", "percent.mito"))
# filteredAL1 <- ScaleData(object = filteredAL1, vars.to.regress = c("nUMI", "percent.mito"))
# filteredAL2 <- ScaleData(object = filteredAL2, vars.to.regress = c("nUMI", "percent.mito"))
# filteredAL3 <- ScaleData(object = filteredAL3, vars.to.regress = c("nUMI", "percent.mito"))
# filteredCPP1 <- ScaleData(object= filteredCPP1, vars.to.regress = c("nUMI", "percent.mito"))
# filteredCPP3 <- ScaleData(object= filteredCPP3, vars.to.regress = c("nUMI", "percent.mito"))
# filteredCPP4 <- ScaleData(object= filteredCPP4, vars.to.regress = c("nUMI", "percent.mito"))
# AL.combined <- ScaleData(object= AL.combined, vars.to.regress = c("nUMI", "percent.mito", "batch"))
# CPN.combined <- ScaleData(object= CPN.combined, vars.to.regress = c("nUMI", "percent.mito", "batch"))
# CPP.combined <- ScaleData(object= CPP.combined, vars.to.regress = c("nUMI", "percent.mito", "batch"))
# 
# save(filteredCPN2, file= "objects/CPN2_seurat_filt_scale")
# save(filteredCPN3, file= "objects/CPN3_seurat_filt_scale")
# save(filteredCPN4, file= "objects/CPN4_seurat_filt_scale")
# save(filteredAL1, file= "objects/AL1_seurat_filt_scale")
# save(filteredAL2, file= "objects/AL2_seurat_filt_scale")
# save(filteredAL3, file= "objects/AL3_seurat_filt_scale")
# save(filteredCPP1, file= "objects/CPP1_seurat_filt_scale")
# save(filteredCPP3, file= "objects/CPP3_seurat_filt_scale")
# save(filteredCPP4, file= "objects/CPP4_seurat_filt_scale")
# save(AL.combined, file= "objects/AL.combined_seurat_filt_merge_scale_batch")
# save(CPN.combined, file= "objects/CPN.combined_seurat_filt_merge_scale_batch")
# save(CPP.combined, file= "objects/CPP.combined_seurat_filt_merge_scale_batch")

ALcellnames <- colnames(AL.combined@data)
AL1cellnames <- colnames(filteredAL1@data)
AL2cellnames <- colnames(filteredAL2@data)
AL3cellnames <- colnames(filteredAL3@data)
CPPcellnames <- colnames(CPP.combined@data)
CPP1cellnames <- colnames(filteredCPP1@data)
CPP3cellnames <- colnames(filteredCPP3@data)
CPP4cellnames <- colnames(filteredCPP4@data)
CPNcellnames <- colnames(CPN.combined@data)
CPN2cellnames <- colnames(filteredCPN2@data)
CPN3cellnames <- colnames(filteredCPN3@data)
CPN4cellnames <- colnames(filteredCPN4@data)
#multi CCA for conditions ####
# 
# ob.list <- list(AL.combined, CPN.combined, CPP.combined)
# genes.use <- c()
# for (i in 1:length(ob.list)) {
#   genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
# }
# genes.use <- names(which(table(genes.use) > 1))
# for (i in 1:length(ob.list)) {
#   genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
# }
# 
# HAM.combined <- RunMultiCCA(object.list = ob.list, genes.use = genes.use, num.ccs = 30)
# pdf("qc/cca/HAM_seurat_filt_merge_scale_batch_cca_dimplot.pdf")
# p1 <- DimPlot(object = HAM.combined, reduction.use = "cca", group.by = "sample", pt.size = 0.2, do.return = TRUE)
# p2 <- VlnPlot(object = HAM.combined, features.plot = "CC1", group.by = "sample", do.return = TRUE)
# plot_grid(p1, p2)
# dev.off()
# 
# save(HAM.combined, file="objects/HAM_agg_merge_cca_unaligned")
# 
# #see how many canonical correlation components you need
# pdf("qc/cca/HAM_seurat_filt_merge_scale_batch_cca_metagenebicorplot.pdf")
# MetageneBicorPlot(HAM.combined, grouping.var = "sample", dims.eval = 1:30)
# dev.off()
# pdf("qc/cca/HAM_seurat_filt_merge_scale_batch_cca_dimplot.pdf")
# p2 <- DimHeatmap(object = HAM.combined, reduction.type = "cca", cells.use = 2000, dim.use = 1:30, do.balanced = TRUE)
# dev.off()
# 
# save(HAM.combined, file="objects/HAM.combined_unaligned")
# 
# #get rid of  outliers
# HAM.combined <- CalcVarExpRatio(object = HAM.combined, reduction.type = "pca", grouping.var = "sample", dims.use = 1:30)
# 
# #We save the original dataset
# HAM.combined.all.save <- HAM.combined
# 
# # Create the corresponding object by removing specific cells
# HAM.combined <- SubsetData(object = HAM.combined, subset.name = "var.ratio.pca", accept.low = 0.5)
# 
# # # Create Seurat object with only the discarded cells (accept.high parameter)
# HAM.combined.discard <- SubsetData(object = HAM.combined.all.save, subset.name = "var.ratio.pca",  accept.high = 0.5)
# # # cat("Median Kept Genes = ",median(HAM.combined.discard@meta.data[, "nGene"]),"\nMedian Discarded Genes = ",median( HAM.combined.discard@meta.data[, "nGene"]))
# 
# save(HAM.combined.all.save, file="objects/HAM_agg_merge_cca_aligned.all.save")
# save(HAM.combined.discard, file="objects/HAM_agg_merge_cca_aligned.discard")
# 
# #alignment
# HAM.combined <- AlignSubspace(object = HAM.combined, reduction.type = "cca", grouping.var = "condition", dims.align =1:23)
# save(HAM.combined, file="objects/HAM_agg_merge_cca_aligned")
# 
# #clustering
# HAM.combined <- FindClusters(object = HAM.combined, reduction.type = "cca.aligned", dims.use = 1:23, save.SNN = TRUE, resolution = 0.6)
# HAM.combined <- RunTSNE(object = HAM.combined, reduction.use = "cca.aligned", dims.use = 1:23, do.fast = TRUE)

# save(HAM.combined, file="objects/HAM_agg_merge_cca_aligned")

# MyColors<-c(brewer.pal(40,"Paired"),"#000000")
## Modified Tsne ####

library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

## Define Colors
MyColors<-c(brewer.pal(20,"Paired"),"#000000")
MyColors[14]<-"#afaeae"
MyColors[15]<-"#ef32c6"

# # # # #clustering ####
# # for (i in c(0.3, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6)){
HAM.combined <- FindClusters(object = HAM.combined, reduction.type = "cca.aligned", dims.use = 1:23, save.SNN = TRUE, resolution = 0.6, force.recalc = T)

HAM.combined <- RunTSNE(object = HAM.combined, check_duplicates= FALSE, reduction.use = "cca.aligned", dims.use = 1:23, do.fast = TRUE)

pdf("tsne/HAM_seurat_filt_merge_scale_batch_cca_tsne_tree.pdf")
CalcAlignmentMetric(HAM.combined, reduction.use = "cca.aligned", dims.use = 1:23, grouping.var =  "condition")
HAM.combined <- BuildClusterTree(HAM.combined, reorder.numeric = TRUE, do.reorder = TRUE, do.plot = TRUE)
# TSNEPlot(object = HAM.combined)
freq_table <- prop.table(x = table(HAM.combined@ident, 
                                   HAM.combined@meta.data[, "condition"]), 
                         margin = 2)
barplot(height = freq_table, legend.text = T, col= MyColors)
TSNEPlot(object = HAM.combined, pt.size = 0.5, cells.use = ALcellnames, 
         # plot.title ="AL",
         group.by = rev("ident"),
         # cols.use=MyColors,
         do.label = T)
# TSNEPlot(object = HAM.combined, pt.size = 0.5, cells.use = AL1cellnames, plot.title ="AL1")
# TSNEPlot(object = HAM.combined, pt.size = 0.5, cells.use = AL2cellnames, plot.title ="AL2")
# TSNEPlot(object = HAM.combined, pt.size = 0.5, cells.use = AL3cellnames, plot.title ="AL3")
TSNEPlot(object = HAM.combined, pt.size = 0.5, cells.use = CPNcellnames, 
         # plot.title ="CPN",
         group.by = rev("ident"),
         # cols.use=MyColors,
         do.label = T)
# TSNEPlot(object = HAM.combined, pt.size = 0.5, cells.use = CPN2cellnames, plot.title ="CPN2")
# TSNEPlot(object = HAM.combined, pt.size = 0.5, cells.use = CPN3cellnames, plot.title ="CPN3")
# TSNEPlot(object = HAM.combined, pt.size = 0.5, cells.use = CPN4cellnames, plot.title ="CPN4")
TSNEPlot(object = HAM.combined, pt.size = 0.5, 
         cells.use = CPPcellnames,
         group.by = rev("ident"),
         plot.title ="CPP",
         # cols.use=MyColors,
         do.label=T)
# TSNEPlot(object = HAM.combined, pt.size = 0.5, cells.use = CPP1cellnames, plot.title ="CPP1")
# TSNEPlot(object = HAM.combined, pt.size = 0.5, cells.use = CPP3cellnames, plot.title ="CPP3")
# TSNEPlot(object = HAM.combined, pt.size = 0.5, cells.use = CPP4cellnames, plot.title ="CPP4")
dev.off()

# plot_grid(TsneCondition, TsneCluster)
# dev.off()
## heatmap ####
# markers <- FindAllMarkers(HAM.combined, only.pos = T, min.pct = 0.25, thresh.use=0.25)
# write.csv(markers, file = "de_tables/HAM_agg_merge_cca_aligned_allmarkers.csv", append = T, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
markers <- read.csv("de_tables/HAM_agg_merge_cca_aligned_allmarkers.csv", header=T, row.names = 1)

MyGenes<-c()
Cluster<-c()
for(jj in unique(markers$cluster))
{
  MK<-markers[markers$cluster==jj,]
  # MK<-MK[order(MK$Specific,decreasing=T),]
  if(nrow(MK)>30)
  {
    MyGenes<-c(MyGenes,as.vector(MK$gene))
    Cluster<-c(Cluster,rep(jj,30))
  }else{

    MyGenes<-c(MyGenes,as.vector(MK$gene))
    Cluster<-c(Cluster,rep(jj,length(MK$gene)))
  }
}

Data<-HAM.combined@scale.data[MyGenes,]
colannotations<-data.frame(Cluster=as.vector(HAM.combined@meta.data$IdentLetter),Condition=as.vector(HAM.combined@meta.data$condition))
rownames(colannotations)<-colnames(Data)
colannotations<-colannotations[order(colannotations$Cluster,colannotations$Condition),,drop=F]
Data<-Data[,match(rownames(colannotations),colnames(Data))]
Data[which(Data>1.5)]<-1.5
Data[which(Data<(-1.5))]<--1.5

hmcols<-colorRampPalette(c("blue","lightgrey","red"))(256)
pheatcolors<-list(Cluster=c("A"=MyColors[1],"B"=MyColors[2],"C"=MyColors[3],"D"=MyColors[4],
                            "E"=MyColors[5],"F"=MyColors[6],"G"=MyColors[7],"H"=MyColors[8],
                            "I"=MyColors[9],"J"=MyColors[10],"K"=MyColors[11], "L"=MyColors[12], "M"=MyColors[13], "N"=MyColors[14], "O"=MyColors[15]),Condition=c("AL"="green","CPN"="blue", "CPP"="red"))
colannotations[,2]<-as.vector(colannotations[,2])
colannotations[,2]<-factor(colannotations[,2])
colannotations<-colannotations[,c(2,1)]
rowannotations<-data.frame(Cluster=Cluster)
rowannotations[,1]<-LETTERS[as.numeric(as.vector(rowannotations[,1]))+1]

rownames(Data)<-paste(MyGenes,"_",Cluster,sep="")
rownames(rowannotations)<-paste(MyGenes,"_",Cluster,sep="")

pdf(file="heatmaps/ham_seurat_monocle_allmarkers.pdf",width=30, height = 20)
# tiff(file="heatmaps/ham_seurat_monocle_allmarkers.tiff",units="in",width=10,height=7,res=150)
pheatmap(Data,annotation_col=colannotations,annotation_row = rowannotations,show_colnames = FALSE,show_rownames = FALSE,cluster_cols = FALSE,
         cluster_rows = FALSE,color =hmcols,annotation_colors = pheatcolors,gaps_col = cumsum(table(colannotations[,2])),
         border_color = "black")
dev.off()

save(HAM.combined, file="objects/HAM_agg_merge_cca_aligned")

# # multi CCA for CPN ####
# 
# ob.list <- list(filteredCPN2, filteredCPN3, filteredCPN4)
# genes.use <- c()
# for (i in 1:length(ob.list)) {
#   genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
# }
# genes.use <- names(which(table(genes.use) > 1))
# for (i in 1:length(ob.list)) {
#   genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
# }
# 
# cpn1234 <- RunMultiCCA(object.list = ob.list, genes.use = genes.use, num.ccs = 30)
# pdf("qc/cca/cpn234_agg_cca.pdf")
# p1 <- DimPlot(object = cpn1234, reduction.use = "cca", group.by = "sample", pt.size = 0.2, do.return = TRUE)
# p2 <- VlnPlot(object = cpn1234, features.plot = "CC1", group.by = "sample", do.return = TRUE)
# plot_grid(p1, p2)
# dev.off()
# 
# save(cpn1234, file="objects/cpn234_agg_cca_unaligned")
# 
# #see how many canonical correlation components you need
# pdf("qc/components_cpn234_agg_cca.pdf")
# MetageneBicorPlot(cpn1234, grouping.var = "sample", dims.eval = 1:30)
# dev.off()
# pdf("qc/heatmapCC_cpn234_agg_cca.pdf")
# p2 <- DimHeatmap(object = cpn1234, reduction.type = "cca", cells.use = 2000, dim.use = 1:30, do.balanced = TRUE)
# dev.off()
# 
# #get rid of  outliers
# cpn1234 <- CalcVarExpRatio(object = cpn1234, reduction.type = "pca", grouping.var = "sample", dims.use = 1:30)
# 
# #We save the original dataset
# cpn1234.all.save <- cpn1234
# 
# # Create the corresponding object by removing specific cells
# cpn1234 <- SubsetData(object = cpn1234, subset.name = "var.ratio.pca", accept.low = 0.5)
# 
# # # Create Seurat object with only the discarded cells (accept.high parameter)
# cpn1234.discard <- SubsetData(object = cpn1234.all.save, subset.name = "var.ratio.pca",  accept.high = 0.5)
# # # cat("Median Kept Genes = ",median(cpn1234.discard@meta.data[, "nGene"]),"\nMedian Discarded Genes = ",median( cpn1234.discard@meta.data[, "nGene"]))
# 
# save(cpn1234.all.save, file="objects/cpn234_agg_cca.all.save")
# save(cpn1234.discard, file="objects/cpn234_agg_cca.discard")
# 
# #alignment
# cpn1234 <- AlignSubspace(object = cpn1234, reduction.type = "cca", grouping.var = "sample", dims.align =1:20)
# cpn1234 <- cpn1234
# save(cpn1234, file="objects/cpn234_cca")
# 
# # #clustering
# cpn1234 <- FindClusters(object = cpn1234, reduction.type = "cca.aligned", dims.use = 1:20, save.SNN = TRUE, resolution = 0.6)
# cpn1234 <- RunTSNE(object = cpn1234, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = TRUE)
# # #
# CPNcellnames <- grep(pattern= "^CPN", x=colnames(cpn1234@data), value=TRUE)
# CPN2cellnames <- grep(pattern= "^CPN2", x=colnames(cpn1234@data), value=TRUE)
# CPN3cellnames <- grep(pattern= "^CPN3", x=colnames(cpn1234@data), value=TRUE)
# CPN4cellnames <- grep(pattern= "^CPN4", x=colnames(cpn1234@data), value=TRUE)
# 
# pdf("tsne/cpn234_agg_cca.pdf")
# TSNEPlot(object = cpn1234, group.by = "sample", do.return = TRUE, pt.size = 0.5)
# TSNEPlot(object = cpn1234, do.return = TRUE, pt.size = 0.5)
# TSNEPlot(object = cpn1234, cells.use = CPN2cellnames, do.return = TRUE, pt.size = 0.5)
# TSNEPlot(object = cpn1234, cells.use = CPN3cellnames, do.return = TRUE, pt.size = 0.5)
# TSNEPlot(object = cpn1234, cells.use = CPN4cellnames, do.return = TRUE, pt.size = 0.5)
# # cpn1234 <- DoKMeans(cpn1234, k.genes = 500, k.cells= 10)
# # TSNEPlot(object = cpn1234, do.return = TRUE, pt.size = 0.5)
# dev.off()
# 
# markers <- FindAllMarkers(cpn1234, only.pos = T, min.pct = 0.25, thresh.use=0.25)
# write.csv(markers, file = "de_tables/CPN_cca_aligned_allmarkers.csv", append = T, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# read.csv("de_tables/CPN_cca_aligned_allmarkers.csv", header=T, row.names = 1)
# 
# # ## Add cluster letters
# cpn1234@meta.data$IdentLetter<-as.factor(LETTERS[as.numeric(cpn1234@meta.data$res.0.6)+1])
# 
# MyGenes<-c()
# Cluster<-c()
# for(jj in unique(markers$cluster))
# {
#   MK<-markers[markers$cluster==jj,]
#   MK<-MK[order(MK$Specific,decreasing=T),]
#   if(nrow(MK)>200)
#   {
#     MyGenes<-c(MyGenes,as.vector(MK$gene))
#     Cluster<-c(Cluster,rep(jj,200))
#   }else{
#     
#     MyGenes<-c(MyGenes,as.vector(MK$gene))
#     Cluster<-c(Cluster,rep(jj,length(MK$gene)))
#   }
# }
# 
# Data<-cpn1234@scale.data[MyGenes,]
# colannotations<-data.frame(Cluster=as.vector(cpn1234@meta.data$IdentLetter),Condition=as.vector(cpn1234@meta.data$condition))
# rownames(colannotations)<-colnames(Data)
# colannotations<-colannotations[order(colannotations$Cluster,colannotations$Condition),,drop=F]
# Data<-Data[,match(rownames(colannotations),colnames(Data))]
# Data[which(Data>1.5)]<-1.5
# Data[which(Data<(-1.5))]<--1.5
# 
# hmcols<-colorRampPalette(c("blue","lightgrey","red"))(256)
# pheatcolors<-list(Cluster=c("A"=MyColors[1],"B"=MyColors[2],"C"=MyColors[3],"D"=MyColors[4],
#                             "E"=MyColors[5],"F"=MyColors[6],"G"=MyColors[7],"H"=MyColors[8],
#                             "I"=MyColors[9],"J"=MyColors[10],"K"=MyColors[11], "L"=MyColors[12], "M"=MyColors[13], "N"=MyColors[14], "O"=MyColors[15]),Condition=c("AL"="green","CPN"="blue", "CPP"="red"))
# colannotations[,2]<-as.vector(colannotations[,2])
# colannotations[,2]<-factor(colannotations[,2])
# colannotations<-colannotations[,c(2,1)]
# rowannotations<-data.frame(Cluster=Cluster)
# rowannotations[,1]<-LETTERS[as.numeric(as.vector(rowannotations[,1]))+1]
# 
# rownames(Data)<-paste(MyGenes,"_",Cluster,sep="")
# rownames(rowannotations)<-paste(MyGenes,"_",Cluster,sep="")
# 
# pdf(file="heatmaps/ham_seurat_monocle_allmarkers.pdf",width=20)
# # tiff(file="heatmaps/ham_seurat_monocle_allmarkers.tiff",units="in",width=10,height=7,res=150)
# pheatmap(Data,annotation_col=colannotations,annotation_row = rowannotations,show_colnames = FALSE,show_rownames = FALSE,cluster_cols = FALSE,
#          cluster_rows = FALSE,color =hmcols,annotation_colors = pheatcolors,gaps_col = cumsum(table(colannotations[,2])),
#          border_color = "black")
# dev.off()
# 
# 
# #multi CCA for CPP
# 
# ob.list <- list(filteredCPP1, filteredCPP3, filteredCPP4)
# genes.use <- c()
# for (i in 1:length(ob.list)) {
#   genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
# }
# genes.use <- names(which(table(genes.use) > 1))
# for (i in 1:length(ob.list)) {
#   genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
# }
# 
# CPP1234 <- RunMultiCCA(object.list = ob.list, genes.use = genes.use, num.ccs = 30)
# pdf("qc/CPP134_agg_cca.pdf")
# p1 <- DimPlot(object = CPP1234, reduction.use = "cca", group.by = "sample", pt.size = 0.2, do.return = TRUE)
# p2 <- VlnPlot(object = CPP1234, features.plot = "CC1", group.by = "sample", do.return = TRUE)
# plot_grid(p1, p2)
# dev.off()
# 
# save(CPP1234, file="objects/CPP134_agg_cca_unaligned")
# 
# #see how many canonical correlation components you need
# pdf("qc/components_CPP134_agg_cca.pdf")
# MetageneBicorPlot(CPP1234, grouping.var = "sample", dims.eval = 1:30)
# dev.off()
# pdf("qc/heatmapCC_CPP134_agg_cca.pdf")
# p2 <- DimHeatmap(object = CPP1234, reduction.type = "cca", cells.use = 2000, dim.use = 1:30, do.balanced = TRUE)
# dev.off()
# 
# #get rid of  outliers
# CPP1234 <- CalcVarExpRatio(object = CPP1234, reduction.type = "pca", grouping.var = "sample", dims.use = 1:30)
# 
# #We save the original dataset
# CPP1234.all.save <- CPP1234
# 
# # Create the corresponding object by removing specific cells
# CPP1234 <- SubsetData(object = CPP1234, subset.name = "var.ratio.pca", accept.low = 0.5)
# 
# # # Create Seurat object with only the discarded cells (accept.high parameter)
# CPP1234.discard <- SubsetData(object = CPP1234.all.save, subset.name = "var.ratio.pca",  accept.high = 0.5)
# # # cat("Median Kept Genes = ",median(CPP1234.discard@meta.data[, "nGene"]),"\nMedian Discarded Genes = ",median( CPP1234.discard@meta.data[, "nGene"]))
# 
# save(CPP1234.all.save, file="objects/CPP134_agg_cca.all.save")
# save(CPP1234.discard, file="objects/CPP134_agg_cca.discard")
# 
# # alignment
# CPP1234 <- AlignSubspace(object = CPP1234, reduction.type = "cca", grouping.var = "sample", dims.align =1:20)
# save(CPP1234, file="objects/CPP134_agg_cca")
# 
# #clustering
# CPP1234 <- FindClusters(object = CPP1234, reduction.type = "cca.aligned", dims.use = 1:20, save.SNN = TRUE, resolution = 0.6)
# CPP1234 <- RunTSNE(object = CPP1234, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = TRUE)
# 
# CPPcellnames <- grep(pattern= "^CPP", x=colnames(CPP1234@data), value=TRUE)
# CPP1cellnames <- grep(pattern= "^CPP1", x=colnames(CPP1234@data), value=TRUE)
# CPP3cellnames <- grep(pattern= "^CPP3", x=colnames(CPP1234@data), value=TRUE)
# CPP4cellnames <- grep(pattern= "^CPP4", x=colnames(CPP1234@data), value=TRUE)
# 
# pdf("tsne/CPP1234_aligned_10dim.pdf")
# TSNEPlot(object = CPP1234, group.by = "sample", do.return = TRUE, pt.size = 0.5)
# TSNEPlot(object = CPP1234, do.return = TRUE, pt.size = 0.5)
# TSNEPlot(object = CPP1234, cells.use = CPP1cellnames, do.return = TRUE, pt.size = 0.5)
# TSNEPlot(object = CPP1234, cells.use = CPP3cellnames, do.return = TRUE, pt.size = 0.5)
# TSNEPlot(object = CPP1234, cells.use = CPP4cellnames, do.return = TRUE, pt.size = 0.5)
# dev.off()
# 
# markers <- FindAllMarkers(CPP1234, only.pos = T, min.pct = 0.25, thresh.use=0.25)
# write.csv(markers, file = "de_tables/CPP_cca_aligned_allmarkers.csv", append = T, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# 
# ## multi CCA for AL
# ob.list <- list(filteredAL1, filteredAL2, filteredAL3)
# genes.use <- c()
# for (i in 1:length(ob.list)) {
#   genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
# }
# genes.use <- names(which(table(genes.use) > 1))
# for (i in 1:length(ob.list)) {
#   genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
# }
# 
# AL123 <- RunMultiCCA(object.list = ob.list, genes.use = genes.use, num.ccs = 30)
# pdf("qc/AL123_agg_cca.pdf")
# p1 <- DimPlot(object = AL123, reduction.use = "cca", group.by = "sample", pt.size = 0.2, do.return = TRUE)
# p2 <- VlnPlot(object = AL123, features.plot = "CC1", group.by = "sample", do.return = TRUE)
# plot_grid(p1, p2)
# dev.off()
# 
# save(AL123, file="objects/AL123_agg_cca_unaligned")
# 
# #see how many canonical correlation components you need
# pdf("qc/components_AL123_agg_cca.pdf")
# MetageneBicorPlot(AL123, grouping.var = "sample", dims.eval = 1:30)
# dev.off()
# pdf("qc/heatmapCC_AL123_agg_cca.pdf")
# p2 <- DimHeatmap(object = AL123, reduction.type = "cca", cells.use = 2000, dim.use = 1:30, do.balanced = TRUE)
# dev.off()
# 
# #get rid of  outliers
# AL123 <- CalcVarExpRatio(object = AL123, reduction.type = "pca", grouping.var = "sample", dims.use = 1:30)
# 
# #We save the original dataset
# AL123.all.save <- AL123
# 
# # Create the corresponding object by removing specific cells
# AL123 <- SubsetData(object = AL123, subset.name = "var.ratio.pca", accept.low = 0.5)
# 
# # # Create Seurat object with only the discarded cells (accept.high parameter)
# AL123.discard <- SubsetData(object = AL123.all.save, subset.name = "var.ratio.pca",  accept.high = 0.5)
# # # cat("Median Kept Genes = ",median(AL123.discard@meta.data[, "nGene"]),"\nMedian Discarded Genes = ",median( AL123.discard@meta.data[, "nGene"]))
# 
# save(AL123.all.save, file="objects/AL123_agg_cca.all.save")
# save(AL123.discard, file="objects/AL123_agg_cca.discard")
# #alignment
# AL123 <- AlignSubspace(object = AL123, reduction.type = "cca", grouping.var = "sample", dims.align =1:20)
# save(AL123, file="objects/AL123_agg_cca")
# 
# clustering
# AL123 <- FindClusters(object = AL123, reduction.type = "cca.aligned", dims.use = 1:20, save.SNN = TRUE, resolution = 0.6)
# AL123 <- RunTSNE(object = AL123, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = TRUE)
# 
# ALcellnames <- grep(pattern= "^AL", x=colnames(AL123@data), value=TRUE)
# AL1cellnames <- grep(pattern= "^AL1", x=colnames(AL123@data), value=TRUE)
# AL2cellnames <- grep(pattern= "^AL2", x=colnames(AL123@data), value=TRUE)
# AL3cellnames <- grep(pattern= "^AL3", x=colnames(AL123@data), value=TRUE)
# 
# pdf("tsne/AL123_aligned_10dim.pdf")
# TSNEPlot(object = AL123, group.by = "sample", do.return = TRUE, pt.size = 0.5)
# TSNEPlot(object = AL123, do.return = TRUE, pt.size = 0.5)
# TSNEPlot(object = AL123, cells.use = AL1cellnames, do.return = TRUE, pt.size = 0.5)
# TSNEPlot(object = AL123, cells.use = AL2cellnames, do.return = TRUE, pt.size = 0.5)
# TSNEPlot(object = AL123, cells.use = AL3cellnames, do.return = TRUE, pt.size = 0.5)
# # AL123 <- DoKMeans(AL123, k.genes = 500, k.cells= 10)
# # TSNEPlot(object = AL123, do.return = TRUE, pt.size = 0.5)
# dev.off()
# 
# markers <- FindAllMarkers(AL123, only.pos = T, min.pct = 0.25, thresh.use=0.25)
# write.csv(markers, file = "de_tables/AL_cca_aligned_allmarkers.csv", append = T, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# 
# ####clustering and tSNE combined elements separated conditions
# # # # AL
# AL.combined <- RunPCA(object = AL.combined, pc.genes = hvg.AL.combined, do.print = T, pcs.print = 1:10, genes.print = 10)
# pdf("qc/AL_agg_merge.pdf")
# VizPCA(object = AL.combined, pcs.use= 1:2)
# PCAPlot(object= AL.combined, dim.1 = 1, dim.2 = 2)
# PCHeatmap(object= AL.combined, pc.use = 1:9, cells.use = 500, do.balanced = T, label.columns = F)
# # AL.combined <- JackStraw(object=AL.combined, num.replicate = 100, do.print=F)
# # JackStrawPlot(object=AL.combined, PCs=7:16)
# PCElbowPlot(object = AL.combined)
# dev.off()
# 
# AL.combined <- FindClusters(object = AL.combined, reduction.type = "pca", save.SNN = TRUE, resolution=0.6, dims.use=1:11)
# AL.combined <- RunTSNE(object = AL.combined, check_duplicates= FALSE, reduction.use = "pca", dims.use = 1:11, do.fast = TRUE)
# pdf("tsne/ALL_agg_merge.pdf")
# TSNEPlot(object = AL.combined)
# TSNEPlot(object = AL.combined, group.by = "sample", colors.use = c("green", "blue", "red"))
# TSNEPlot(object = AL.combined, pt.size = 0.5, cells.use = AL1cellnames, plot.title ="AL1")
# TSNEPlot(object = AL.combined, pt.size = 0.5, cells.use = AL2cellnames, plot.title ="AL2")
# TSNEPlot(object = AL.combined, pt.size = 0.5, cells.use = AL3cellnames, plot.title ="AL3")
# dev.off()
# 
# save(AL.combined, file= "objects/ALL_agg_merge")
# 
# markers <- FindAllMarkers(AL.combined, only.pos = T, min.pct = 0.25, thresh.use=0.25)
# write.csv(markers, file = "de_tables/AL_merge_batch_allmarkers.csv", append = T, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# 
# # # # CPN
# CPN.combined <- RunPCA(object = CPN.combined, pc.genes = hvg.CPN.combined, do.print = T, pcs.print = 1:10, genes.print = 10)
# pdf("qc/CPNL_agg_merge.pdf")
# VizPCA(object = CPN.combined, pcs.use= 1:2)
# PCAPlot(object= CPN.combined, dim.1 = 1, dim.2 = 2)
# PCHeatmap(object= CPN.combined, pc.use = 1:9, cells.use = 500, do.balanced = T, label.columns = F)
# # CPN.combined <- JackStraw(object=CPN.combined, num.replicate = 100, do.print=F)
# # JackStrawPlot(object=CPN.combined, PCs=7:16)
# PCElbowPlot(object = CPN.combined)
# dev.off()
# 
# CPN.combined <- FindClusters(object = CPN.combined, reduction.type = "pca", save.SNN = TRUE, resolution=0.6, dims.use=1:11)
# CPN.combined <- RunTSNE(object = CPN.combined, check_duplicates= FALSE, reduction.use = "pca", dims.use = 1:11, do.fast = TRUE)
# pdf("tsne/CPNL_agg_merge.pdf")
# TSNEPlot(object = CPN.combined)
# TSNEPlot(object = CPN.combined, group.by = "sample", colors.use = c("green", "blue", "red"))
# TSNEPlot(object = CPN.combined, pt.size = 0.5, cells.use = CPN2cellnames, plot.title ="CPN2")
# TSNEPlot(object = CPN.combined, pt.size = 0.5, cells.use = CPN3cellnames, plot.title ="CPN3")
# TSNEPlot(object = CPN.combined, pt.size = 0.5, cells.use = CPN4cellnames, plot.title ="CPN4")
# dev.off()
# 
# save(CPN.combined, file= "objects/CPN_agg_merge")
# 
# markers <- FindAllMarkers(CPN.combined, only.pos = T, min.pct = 0.25, thresh.use=0.25)
# write.csv(markers, file = "de_tables/CPN_merge_batch_allmarkers.csv", append = T, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# 
# # # # CPP
# CPP.combined <- RunPCA(object = CPP.combined, pc.genes = hvg.CPP.combined, do.print = T, pcs.print = 1:10, genes.print = 10)
# pdf("qc/CPP_agg_merge.pdf")
# VizPCA(object = CPP.combined, pcs.use= 1:2)
# PCAPlot(object= CPP.combined, dim.1 = 1, dim.2 = 2)
# PCHeatmap(object= CPP.combined, pc.use = 1:9, cells.use = 500, do.balanced = T, label.columns = F)
# # CPP.combined <- JackStraw(object=CPP.combined, num.replicate = 100, do.print=F)
# # JackStrawPlot(object=CPP.combined, PCs=7:16)
# PCElbowPlot(object = CPP.combined)
# dev.off()
# 
# CPP.combined <- FindClusters(object = CPP.combined, reduction.type = "pca", save.SNN = TRUE, resolution=0.6, dims.use=1:11)
# CPP.combined <- RunTSNE(object = CPP.combined, check_duplicates= FALSE, reduction.use = "pca", dims.use = 1:11, do.fast = TRUE)
# pdf("tsne/CPP_agg_merge.pdf")
# TSNEPlot(object = CPP.combined)
# TSNEPlot(object = CPP.combined, group.by = "sample", colors.use = c("green", "blue", "red"))
# TSNEPlot(object = CPP.combined, pt.size = 0.5, cells.use = CPP1cellnames, plot.title ="CPP1")
# TSNEPlot(object = CPP.combined, pt.size = 0.5, cells.use = CPP3cellnames, plot.title ="CPP3")
# TSNEPlot(object = CPP.combined, pt.size = 0.5, cells.use = CPP4cellnames, plot.title ="CPP4")
# dev.off()
# 
# save(CPP.combined, file= "objects/CPP_agg_merge")
# 
# markers <- FindAllMarkers(CPP.combined, only.pos = T, min.pct = 0.25, thresh.use=0.25)
# write.csv(markers, file = "de_tables/CPP_merge_batch_allmarkers.csv", append = T, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# 
# # MyColors<-c(brewer.pal(40,"Paired"),"#000000")
# ## Modified Tsne ####
# 
# library(pheatmap)
# library(RColorBrewer)
# library(ggplot2)
# library(dplyr)
# 
# ## Define Colors
# MyColors<-c(brewer.pal(20,"Paired"),"#000000")
# MyColors[8]<-"#afaeae"
# MyColors[11]<-"#ef32c6"
# 
# # # # #clustering ####
# pdf("qc/dotplot_mergeHAM_0.6.pdf")
# # for (i in c(0.3, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6)){
# cpn1234 <- FindClusters(object = cpn1234, reduction.type = "cca.aligned", dims.use = 1:20, save.SNN = TRUE, resolution = 0.6, force.recalc = T)
# freq_table <- prop.table(x = table(cpn1234@ident, cpn1234@meta.data[, "condition"]),
#                          margin = 2)
# print(freq_table)
# barplot(height = freq_table, col=MyColors, legend.text = T)
# # }
# # dev.off()
# 
# cpn1234 <- RunTSNE(object = cpn1234, check_duplicates= FALSE, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = TRUE)
# 
# pdf("tsne/ham_agg_cca_tree.pdf")
# CalcAlignmentMetric(cpn1234, reduction.use = "cca.aligned", dims.use = 1:20, grouping.var =  "condition")
# cpn1234 <- BuildClusterTree(cpn1234, reorder.numeric = TRUE, do.reorder = TRUE, do.plot = TRUE)
# TSNEPlot(object = cpn1234)
# # TSNEPlot(object = cpn1234, group.by = "condition", colors.use = c("green", "blue", "red"))
# TSNEPlot(object = cpn1234, pt.size = 0.5, cells.use = ALcellnames, plot.title ="AL")
# TSNEPlot(object = cpn1234, pt.size = 0.5, cells.use = AL1cellnames, plot.title ="AL1")
# TSNEPlot(object = cpn1234, pt.size = 0.5, cells.use = AL2cellnames, plot.title ="AL2")
# TSNEPlot(object = cpn1234, pt.size = 0.5, cells.use = AL3cellnames, plot.title ="AL3")
# TSNEPlot(object = cpn1234, pt.size = 0.5, cells.use = CPNcellnames, plot.title ="CPN")
# TSNEPlot(object = cpn1234, pt.size = 0.5, cells.use = CPN2cellnames, plot.title ="CPN2")
# TSNEPlot(object = cpn1234, pt.size = 0.5, cells.use = CPN3cellnames, plot.title ="CPN3")
# TSNEPlot(object = cpn1234, pt.size = 0.5, cells.use = CPN4cellnames, plot.title ="CPN4")
# TSNEPlot(object = cpn1234, pt.size = 0.5, cells.use = CPPcellnames, plot.title ="CPP")
# TSNEPlot(object = cpn1234, pt.size = 0.5, cells.use = CPP1cellnames, plot.title ="CPP1")
# TSNEPlot(object = cpn1234, pt.size = 0.5, cells.use = CPP3cellnames, plot.title ="CPP3")
# TSNEPlot(object = cpn1234, pt.size = 0.5, cells.use = CPP4cellnames, plot.title ="CPP4")
# dev.off()
# 
# # Modified Tsne ####
# 
# library(pheatmap)
# library(RColorBrewer)
# library(ggplot2)
# library(dplyr)
# 
# ## Define Colors
# MyColors<-c(brewer.pal(20,"Paired"),"#000000")
# MyColors[8]<-"#afaeae"
# MyColors[11]<-"#ef32c6"
# 
# ## Add cluster letters
# HAM.combined@meta.data$IdentLetter<-as.factor(LETTERS[as.numeric(HAM.combined@meta.data$res.0.4)+1])
# 
# # MarkersAlign<-cbind(MarkersAlign,MarkersAlign$pct.1-MarkersAlign$pct.2)
# # colnames(MarkersAlign)[8]<-"Specific"
# 
# ## Modified Tsne
# 
# dim(HAM.combined@dr$tsne@cell.embeddings)
# TsneData<-data.frame(x=HAM.combined@dr$tsne@cell.embeddings[,1],y=HAM.combined@dr$tsne@cell.embeddings[,2],Condition=HAM.combined@meta.data$condition,Cluster=HAM.combined@meta.data$IdentLetter)
# Cond<-as.vector(TsneData[,3])
# Cond<-as.factor(Cond)
# TsneData[,3]<-Cond
# 
# centers <- TsneData %>% dplyr::group_by(Cluster) %>% summarize(x = median(x = x),y = median(x = y))
# 
# pdf("tsne/ham_seurat_merge_cca")
# TsneCondition <- ggplot(data = TsneData, mapping = aes(x = x, y = y)) +
#   geom_point(mapping = aes(colour = factor(x = Condition)),size = 0.3,alpha=0.6) +
#   scale_colour_manual(values = c("#ff0b00","#00c4ff", "#025802")) +
#   xlab(label = "t-SNE 1") + ylab(label = "t-SNE 2") + xlim(-40, 40) + ylim(-40, 40)+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank())+
#   guides(colour = guide_legend(override.aes = list(size = 3)))
# 
# 
# TsneCluster<- ggplot(data = TsneData, mapping = aes(x = x, y = y)) +
#   geom_point(mapping = aes(colour = factor(x = Cluster)),size = 0.3) +
#   scale_colour_manual(values = MyColors) +
#   xlab(label = "t-SNE 1") + ylab(label = "t-SNE 2") + xlim(-40, 40) + ylim(-45, 45)+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank())+
#   guides(colour = guide_legend(override.aes = list(size = 3))) +
#   geom_point(data = centers, mapping = aes(x = x, y = y), size = 0, alpha = 0) +
#   geom_text(data = centers,  mapping = aes(label = Cluster), size = 4,fontface = "bold")
# 
# plot_grid(TsneCondition, TsneCluster)
# dev.off()
# ## heatmap ####
# markers <- FindAllMarkers(HAM.combined, only.pos = T, min.pct = 0.25, thresh.use=0.25)
# write.csv(markers, file = "de_tables/mergeHAM_seurat_merge_cca_allmarkers.csv", append = T, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# 
# # ## Add cluster letters
# HAM.combined@meta.data$IdentLetter<-as.factor(LETTERS[as.numeric(HAM.combined@meta.data$res.0.6)+1])
# 
# MyGenes<-c()
# Cluster<-c()
# for(jj in unique(markers$cluster))
# {
#   MK<-markers[markers$cluster==jj,]
#   # MK<-MK[order(MK$Specific,decreasing=T),]
#   if(nrow(MK)>30)
#   {
#     MyGenes<-c(MyGenes,as.vector(MK$gene))
#     Cluster<-c(Cluster,rep(jj,30))
#   }else{
#     
#     MyGenes<-c(MyGenes,as.vector(MK$gene))
#     Cluster<-c(Cluster,rep(jj,length(MK$gene)))
#   }
# }
# 
# Data<-HAM.combined@scale.data[MyGenes,]
# colannotations<-data.frame(Cluster=as.vector(HAM.combined@meta.data$IdentLetter),Condition=as.vector(HAM.combined@meta.data$condition))
# rownames(colannotations)<-colnames(Data)
# colannotations<-colannotations[order(colannotations$Cluster,colannotations$Condition),,drop=F]
# Data<-Data[,match(rownames(colannotations),colnames(Data))]
# Data[which(Data>1.5)]<-1.5
# Data[which(Data<(-1.5))]<--1.5
# 
# hmcols<-colorRampPalette(c("blue","lightgrey","red"))(256)
# pheatcolors<-list(Cluster=c("A"=MyColors[1],"B"=MyColors[2],"C"=MyColors[3],"D"=MyColors[4],
#                             "E"=MyColors[5],"F"=MyColors[6],"G"=MyColors[7],"H"=MyColors[8],
#                             "I"=MyColors[9],"J"=MyColors[10],"K"=MyColors[11], "L"=MyColors[12], "M"=MyColors[13]),Condition=c("AL"="green","CPN"="blue", "CPP"="red"))
# colannotations[,2]<-as.vector(colannotations[,2])
# colannotations[,2]<-factor(colannotations[,2])
# colannotations<-colannotations[,c(2,1)]
# rowannotations<-data.frame(Cluster=Cluster)
# rowannotations[,1]<-LETTERS[as.numeric(as.vector(rowannotations[,1]))+1]
# 
# rownames(Data)<-paste(MyGenes,"_",Cluster,sep="")
# rownames(rowannotations)<-paste(MyGenes,"_",Cluster,sep="")
# 
# pdf(file="heatmaps/ham_seurat_monocle_allmarkers.pdf",width=20)
# # tiff(file="heatmaps/ham_seurat_monocle_allmarkers.tiff",units="in",width=10,height=7,res=150)
# pheatmap(Data,annotation_col=colannotations,annotation_row = rowannotations,show_colnames = FALSE,show_rownames = FALSE,cluster_cols = FALSE,
#          cluster_rows = FALSE,color =hmcols,annotation_colors = pheatcolors,gaps_col = cumsum(table(colannotations[,2])),
#          border_color = "black")
# dev.off()
# 
# one.markers <- FindConservedMarkers(HAM.combined, ident.1 = 1, grouping.var = "condition", print.bar = T)
# write.csv(one.markers, file = "de_tables/HAM.combined2_seurat_merge_cca_consmark_cl1.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# two.markers <- FindConservedMarkers(HAM.combined, ident.1 = 2, grouping.var = "condition", print.bar = T)
# write.csv(two.markers, file = "de_tables/HAM.combined2_seurat_merge_cca_consmark_cl2.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# three.markers <- FindConservedMarkers(HAM.combined, ident.1 = 3, grouping.var = "condition", print.bar = T)
# write.csv(three.markers, file = "de_tables/HAM.combined2_seurat_merge_cca_consmark_cl3.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# four.markers <- FindConservedMarkers(HAM.combined, ident.1 = 4, grouping.var = "condition", print.bar = T)
# write.csv(four.markers, file = "de_tables/HAM.combined2_seurat_merge_cca_consmark_cl4.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# five.markers <- FindConservedMarkers(HAM.combined, ident.1 = 5, grouping.var = "condition", print.bar = T)
# write.csv(five.markers, file = "de_tables/HAM.combined2_seurat_merge_cca_consmark_cl5.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# six.markers <- FindConservedMarkers(HAM.combined, ident.1 = 6, grouping.var = "condition", print.bar = T)
# write.csv(six.markers, file = "de_tables/HAM.combined2_seurat_merge_cca_consmark_cl6.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# seven.markers <- FindConservedMarkers(HAM.combined, ident.1 = 7, grouping.var = "condition", print.bar = T)
# write.csv(seven.markers, file = "de_tables/HAM.combined2_seurat_merge_cca_consmark_cl7.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names= TRUE)
# 
# genes.use = union_all(rownames(one.markers)[1:5], 
#                       rownames(two.markers)[1:5], 
#                       rownames(three.markers)[1:5], 
#                       rownames(four.markers)[1:5], 
#                       rownames(five.markers)[1:5],
#                       rownames(six.markers)[1:5],
#                       rownames(seven.markers)[1:5])
# 
# pdf("tsne/splitdotplot.pdf")
# SplitDotPlotGG(HAM.combined, genes.plot = rev(markers.to.plot), cols.use = c("blue", "red", "green"), x.lab.rot = T, plot.legend = T, dot.scale = 8, do.return = T, grouping.var = "condition")
# dev.off()
# 
# AssessNodes(HAM.combined)
# AssessSplit(HAM.combined, 1, cluster1= 1, cluster2= 2, print.output = TRUE, verbose=T, show.progress = TRUE)
# head(AverageDetectionRate(object = HAM.combined))
# 
# 
# # # # #save objects
# save(HAM.combined, file="objects/seurat_no_agg/mergeHAM")
