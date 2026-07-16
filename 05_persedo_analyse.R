# # 安装Monocle 3
# install.packages("remotes")
# remotes::install_github("cole-trapnell-lab/monocle3")
# install.packages("R.utils")
# remotes::install_github("satijalab/seurat-wrappers")
# 加载库
library(monocle3)
library(SeuratWrappers)  # 用于Seurat和Monocle的结合

#加载数据

load("./subannoComplete.Rdata")



# 假设 subsce 是你的 Seurat 对象
# 将 Seurat 对象转换为 Monocle 3 的 cell_data_set 对象
cds <- as.cell_data_set(subsce)

# 确保亚群注释信息（subcelltype）已包含在 Monocle 对象中
cds@colData$subcelltype <- subsce@meta.data$subcelltype

# 创建细胞聚类和分区
cds <- cluster_cells(cds,reduction_method = 'UMAP')

# 进行 UMAP 降维
cds <- reduce_dimension(cds)

# 构建轨迹图
cds <- learn_graph(cds)

# 找到 Pro-NPCs 细胞的 ID
pro_cells <- colnames(cds)[cds@colData$subcelltype == "pro"]

# 以 Pro-NPCs 细胞为起点进行伪时序排序
cds <- order_cells(cds, root_cells = pro_cells)


# 根据伪时序绘制轨迹图
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE, 
           label_leaves = TRUE, label_branch_points = TRUE)

# 根据 subcelltype 显示轨迹图
plot_cells(cds, color_cells_by = "subcelltype", label_groups_by_cluster = TRUE, 
           label_leaves = TRUE, label_branch_points = TRUE)

# 识别动态变化的基因
deg_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

# 筛选显著差异的基因
deg_genes_sig <- subset(deg_genes, q_value < 0.05)

# 查看前几个动态基因
head(deg_genes_sig)

# 保存伪时序轨迹中的细胞信息
write.csv(as.data.frame(colData(cds)), file = "pseudotime_results.csv", row.names = TRUE)