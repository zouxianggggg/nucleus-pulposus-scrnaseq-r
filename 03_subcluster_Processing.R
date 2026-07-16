library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
load("./NPCells.Rdata")
#对NP细胞再进行亚群注释
#标记基因
np_markers <- list(
  "Pro-NPCs" = c("ACAN", "SOX9", "UBE2C", "TOP2A"),
  "Fibro-NPCs" = c("ACAN", "SOX9", "FBLN1", "FBLN3", "FBLN5"),
  "IR-NPCs" = c("ACAN", "SOX9", "CH3L2"),
  "Met-NPCs" = c("ACAN", "SOX9", "DKK1"),
  "Adh-NPCs" = c("ACAN", "SOX9", "MSMO1"),
  "SR-NPCs" = c("ACAN", "SOX9", "CP")
)

#之前的流程先来一套
scRNAsub <- NormalizeData(scRNAsub)

scRNAsub <- FindVariableFeatures(scRNAsub,selection.method = "vst", nfeatures = 2000)

scale.genes <- rownames(scRNAsub)

scRNAsub <- ScaleData(scRNAsub,features = scale.genes)

scRNAsub <- RunPCA(scRNAsub,features = VariableFeatures(scRNAsub))

DimPlot(scRNAsub, reduction = "pca", group.by = "group")
DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")

ElbowPlot(scRNAsub,ndims = 30,reduction = "pca")

scRNAsub <- FindNeighbors(scRNAsub,dims = 1:20)

scRNAsub <- FindClusters(scRNAsub,resolution = seq(from = 0.6, to = 2.0, by = 0.2))

scRNAsub <- RunUMAP(scRNAsub,dims = 1:20)
scRNAsub <- RunTSNE(scRNAsub,dims = 1:20)

#保存降维聚类数据
save(scRNAsub,file = "./scRNAsub_cluster.Rdata")

load("./scRNAsub_cluster.Rdata")

library(clustree)
#查看聚类树
clustree(scRNAsub)
#Idents(scRNAsub)
##直接选择默认的resolution 1.0

#选择0。4分辨率的cluster
scRNAsub$RNA_snn_res.0.5

Idents(scRNAsub) <- "RNA_snn_res.0.5"

DimPlot(scRNAsub,reduction = "umap", group.by = "orig.ident")


#绘图
umap_integrated1 <- DimPlot(scRNAsub,reduction = "umap", group.by = "orig.ident")
umap_integrated2 <- DimPlot(scRNAsub,reduction = "umap", group.by = "group")
tsne_integrated3 <- DimPlot(scRNAsub,reduction = "tsne", label = TRUE)
umap_integrated4 <- DimPlot(scRNAsub,reduction = "umap", label = TRUE)

umap_tsne_integrated <- CombinePlots(list(umap_integrated1,umap_integrated2,tsne_integrated3,umap_integrated4))

umap_tsne_integrated

#已经是scRNAsub_cluster.Rdata了


