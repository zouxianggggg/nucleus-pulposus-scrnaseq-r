# 加载必要的 R 包
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(harmony)  # Harmony 包
library(tidyverse)
library(patchwork)
library(cowplot)
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
if(!require("BiocManager")) install.packages("BiocManager")
#BiocManager::install('glmGamPoi')
library(glmGamPoi)

# 1. 创建 Seurat 对象
SHIVDD_1 <- CreateSeuratObject(Read10X(data.dir = "../分析/ABFC20240729-01/ABFC20240729-01-SCS-result/2.Cellranger/SHIVDD-1-10XSC3/filtered_feature_bc_matrix"), project = "SHIVDD-1")
SHIVDD_2 <- CreateSeuratObject(Read10X(data.dir = "../分析/ABFC20240729-01/ABFC20240729-01-SCS-result/2.Cellranger/SHIVDD-2-10XSC3/filtered_feature_bc_matrix"), project = "SHIVDD-2")
SHIVDD_3 <- CreateSeuratObject(Read10X(data.dir = "../分析/ABFC20240729-03/ABFC20240729-03-SCS-result/2.Cellranger/SHIVDD-3-10XSC3/filtered_feature_bc_matrix"), project = "SHIVDD-3")
SHCON_1 <- CreateSeuratObject(Read10X(data.dir = "../分析/ABFC20240729-04/ABFC20240729-04-SCS-result/2.Cellranger/SHcon-2-10XSC3/filtered_feature_bc_matrix"), project = "SHCON-1")
SHCON_2 <- CreateSeuratObject(Read10X(data.dir = "../分析/ABFC20240729-06/ABFC20240729-06-SCS-result/2.Cellranger/SHcon-3-10XSC3/filtered_feature_bc_matrix"), project = "SHCON-2")
SHCON_3 <- CreateSeuratObject(Read10X(data.dir = "../分析/ABFC20240729-06/ABFC20240729-06-SCS-result/2.Cellranger/SHcon-4-10XSC3/filtered_feature_bc_matrix"), project = "SHCON-3")

# 2. 计算线粒体基因比例
# 在Seurat V5中，计算线粒体基因比例的函数仍然是 PercentageFeatureSet
# 确保线粒体基因的命名符合 "^MT-" 的模式
SHIVDD_1[["percent.mt"]] <- PercentageFeatureSet(SHIVDD_1, pattern = "^MT-")
SHIVDD_2[["percent.mt"]] <- PercentageFeatureSet(SHIVDD_2, pattern = "^MT-")
SHIVDD_3[["percent.mt"]] <- PercentageFeatureSet(SHIVDD_3, pattern = "^MT-")
SHCON_1[["percent.mt"]] <- PercentageFeatureSet(SHCON_1, pattern = "^MT-")
SHCON_2[["percent.mt"]] <- PercentageFeatureSet(SHCON_2, pattern = "^MT-")
SHCON_3[["percent.mt"]] <- PercentageFeatureSet(SHCON_3, pattern = "^MT-")

# 3. 合并所有样本的 Seurat 对象
combined <- merge(SHIVDD_1, y = c(SHIVDD_2, SHIVDD_3, SHCON_1, SHCON_2, SHCON_3),
                  add.cell.ids = c("SHIVDD-1", "SHIVDD-2", "SHIVDD-3", "SHCON-1", "SHCON-2", "SHCON-3"),
                  project = "Combined")

# 4. 绘制总体的线粒体基因比例和基因数的散点图
FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept = 25, color = "red") + 
  ggtitle("Percent.mt vs nCount_RNA - Combined")

FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = 200, color = "red") +
  ggtitle("nFeature_RNA vs nCount_RNA - Combined")

# 5. 绘制小提琴图，展示基因数、UMI计数、线粒体基因比例的分布
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by = "orig.ident", layer = "count",ncol = 3)

#统计一下细胞数
table(combined[[]]$orig.ident)

# 识别不同组的细胞
# SDD组
SDD_cells <- WhichCells(combined, idents = c("SHIVDD-1", "SHIVDD-2", "SHIVDD-3"))

# MDD组 (剩余的对照组)
MDD_cells <- WhichCells(combined, idents = c("SHCON-1", "SHCON-2", "SHCON-3"))

# 创建 'group' 列，并将分组信息加入到 meta.data 中
combined$group <- NA  # 先创建一个空列

# 给 SDD 组赋值为 "SDD"
combined$group[SDD_cells] <- "SDD"

# 给 MDD 组赋值为 "MDD"
combined$group[MDD_cells] <- "MDD"

# 检查meta.data
head(combined@meta.data)

#归一化数据
combined <- NormalizeData(combined)
#找到高变基因
combined <- FindVariableFeatures(combined)
#回归线粒体基因比例并进行数据缩放
combined <- ScaleData(combined, vars.to.regress = "percent.mt")
#运行 PCA
combined <- RunPCA(combined, verbose = F)

# 可视化 PCA 结果
DimPlot(combined, reduction = "pca", group.by = "group")
DimPlot(combined, reduction = "pca", group.by = "orig.ident")

# 查看高变基因
top10_var_genes <- head(VariableFeatures(combined), 10)
print(top10_var_genes)

#harmony降维
scRNA_harmony <- IntegrateLayers(object = combined,
                                 method = HarmonyIntegration,
                                 orig.reduction = "pca",
                                 new.reduction = "harmony",
                                 verbose = FALSE)
#joinlayers合并样本的data和counts

scRNA_harmony[["RNA"]] <- JoinLayers(scRNA_harmony[["RNA"]])

#降维聚类
#reduction选择"harmony"
ElbowPlot(scRNA_harmony,ndims = 50)
scRNA_harmony <- FindNeighbors(scRNA_harmony,reduction = "harmony",dims = 1:20)

scRNA_harmony <- FindClusters(scRNA_harmony,resolution = seq(from = 0.1, to = 1.0, by = 0.1))

scRNA_harmony <- RunUMAP(scRNA_harmony,dims = 1:20,reduction = "harmony")
scRNA_harmony <- RunTSNE(scRNA_harmony,dims = 1:20,reduction = "harmony")
#保存降维聚类数据
save(scRNA_harmony,file = "scRNA_harmony.Rdata")

library(clustree)
# install.packages("systemfonts", dependencies = TRUE)
# install.packages("igraph")

#查看聚类树
clustree(scRNA_harmony)

#Umap图看下Harmony整合情况
DimPlot(scRNA_harmony,reduction = "umap", group.by = "orig.ident") + ggtitle("Harmony")

#选择0。4分辨率的cluster
scRNA_harmony$RNA_snn_res.0.4

Idents(scRNA_harmony) <- "RNA_snn_res.0.4"

DimPlot(scRNA_harmony,reduction = "umap")

#绘图
umap_integrated1 <- DimPlot(scRNA_harmony,reduction = "umap", group.by = "orig.ident")
umap_integrated2 <- DimPlot(scRNA_harmony,reduction = "umap", group.by = "group")
tsne_integrated3 <- DimPlot(scRNA_harmony,reduction = "tsne", label = TRUE)
umap_integrated4 <- DimPlot(scRNA_harmony,reduction = "umap", label = TRUE)

umap_tsne_integrated <- CombinePlots(list(umap_integrated1,umap_integrated2,tsne_integrated3,umap_integrated4))

umap_tsne_integrated

#MDD SDD对比图
DimPlot(scRNA_harmony,reduction = "umap",split.by = "group")

#保存数据
save(scRNA_harmony,file = "scRNA_harmony2.Rdata")

#查看细胞周期的影响
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scRNA_harmony <- CellCycleScoring(scRNA_harmony,s.features = s.genes,g2m.features = g2m.genes,set.ident = TRUE)
DimPlot(scRNA_harmony,group.by = "Phase")
DimPlot(scRNA_harmony,group.by = "Phase",reduction = "pca")


#提取counts，data数据
#
#提取counts
scRNA_counts = LayerData(scRNA_harmony,assay = "RNA",layer = "counts")
scRNA_counts = as.data.frame(scRNA_counts)
#提取data
scRNA_data = LayerData(scRNA_harmony,assay = "RNA",layer = "data")
scRNA_data = as.data.frame(scRNA_data)

#求每个cluster的平均表达量
av_seurat = AverageExpression(scRNA_harmony,assays = "RNA",group.by = "RNA_snn_res.0.4")
av_seurat = as.data.frame(av_seurat$RNA)
av_seurat
