library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
install.packages('devtools')
devtools::install_github('immunogenomics/presto')
load("./scRNA_harmony2.Rdata")


install.packages("tensorflow")
tensorflow::install_tensorflow(extra_packages='tensorflow-probability')


plot1 = DimPlot(scRNA_harmony,reduction = "umap",label = T)
plot1
plot2 = DimPlot(scRNA_harmony,reduction = "umap",group.by = "orig.ident")
plot2

#tsne
plot3 = DimPlot(scRNA_harmony,reduction = "tsne",label = T)
plot3
plot4 = DimPlot(scRNA_harmony,reduction = "tsne",group.by = "orig.ident")
plot4



#创建一个注释表格
celltype = data.frame(ClusterID = 0:11,
                      celltype = 'unknown')
celltype[celltype$ClusterID %in% c(0,1,2,3,4,5,6,10),2] = "NP cell"
celltype[celltype$ClusterID %in% c(7),2] = "macrophage"
celltype[celltype$ClusterID %in% c(8),2] = "SMC"
celltype[celltype$ClusterID %in% c(9),2] = "Endothelia"
celltype[celltype$ClusterID %in% c(11),2] = "Erythrocyte"

table(celltype$celltype)


#把数据放到另一个对象里
sce.in = scRNA_harmony
sce.in@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.in@meta.data[which(sce.in@meta.data$RNA_snn_res.0.4 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}

#查看每个聚类里每种细胞多少个
table(sce.in@meta.data$celltype,sce.in@meta.data$RNA_snn_res.0.4)
table(sce.in@meta.data$group)

sce = sce.in
DimPlot(sce,reduction = "umap",group.by = "celltype",label = T)
DimPlot(sce,reduction = "tsne",group.by = "celltype",label = T)

#看某个基因的表达情况
VlnPlot(sce,features = "ACAN",group.by = "celltype")

#筛选出NP细胞 
DefaultAssay(sce)
Cells.sub <- subset(sce@meta.data,celltype == "NP cell")
scRNAsub <- subset(sce,cells = row.names(Cells.sub))
#先保存一下，NP细胞集
save(scRNAsub,file = "./NPCells.Rdata")
