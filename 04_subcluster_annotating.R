library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
load("./scRNAsub_cluster.Rdata")

scRNAsub$RNA_snn_res.1.2

Idents(scRNAsub) <- "RNA_snn_res.1.4"

DimPlot(scRNAsub,reduction = "umap", group.by = "orig.ident")
plot1 = DimPlot(scRNAsub,reduction = "umap",label = T)
plot1
plot2 = DimPlot(scRNAsub,reduction = "umap",group.by = "orig.ident")
plot2

#tsne
plot3 = DimPlot(scRNAsub,reduction = "tsne",label = T)
plot3
plot4 = DimPlot(scRNAsub,reduction = "tsne",group.by = "orig.ident")
plot4


np_markers <- list(
  "Pro-NPCs" = c("UBE2C", "TOP2A"),
  "Fibro-NPCs" = c("FBLN1","COL1A1","COL1A2"),
  "IR-NPCs" = c("CHI3L2"),
  "Met-NPCs" = c("DKK1"),
  "Adh-NPCs" = c("MSMO1"),
  "SR-NPCs" = c("CP")
)

gene_tocheck = list("Pro-NPCs" = c("ACAN", "SOX9", "UBE2C", "TOP2A"))
gene_tocheck = list("Fibro-NPCs" = c("ACAN", "SOX9", "FBLN1"))
gene_tocheck = list("Met-NPCs" = c("ACAN", "SOX9", "DKK1"))
gene_tocheck = list("Adh-NPCs" = c("ACAN", "SOX9", "MSMO1"))
gene_tocheck = list("SR-NPCs" = c("ACAN", "SOX9", "CP"))

# 绘制dotplot
DotPlot(scRNAsub, features = np_markers,assay = "RNA",group.by = "RNA_snn_res.1.4") + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




#创建一个注释表格
celltype = data.frame(ClusterID = 0:24,
                      celltype = 'unknown')
celltype[celltype$ClusterID %in% c(12),2] = "pro"
celltype[celltype$ClusterID %in% c(21,16,18,11),2] = "fibro"
#celltype[celltype$ClusterID %in% c(5,13),2] = "met"
#celltype[celltype$ClusterID %in% c(0,3,7,12),2] = "adh"
#celltype[celltype$ClusterID %in% c(1,2,4,8,10),2] = "sr"

table(celltype$celltype)

#把数据放到另一个对象里
subsce.in = scRNAsub
subsce.in@meta.data$subcelltype = "NA"
for(i in 1:nrow(celltype)){
  subsce.in@meta.data[which(subsce.in@meta.data$RNA_snn_res.1.4 == celltype$ClusterID[i]),'subcelltype'] <- celltype$celltype[i]
}

#查看每个聚类里每种细胞多少个
table(subsce.in@meta.data$subcelltype,subsce.in@meta.data$RNA_snn_res.1.4)
table(subsce.in@meta.data$group)

subsce = subsce.in
DimPlot(subsce,reduction = "umap",group.by = "subcelltype",label = T)
DimPlot(subsce,reduction = "tsne",group.by = "subcelltype",label = T)

#保存为新的文件
save(subsce,file = "./subannoComplete.Rdata")
