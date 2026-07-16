library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)

{
  load("./SDDvsMDD差异基因.Rdata")
  load("./fibrovsother差异基因.Rdata")
  
  
  # 确保 markers 数据框按 avg_log2FC 升序排序
  markers_sorted <- top_left %>%
    arrange(avg_log2FC)
  
  # 查找 HDAC9 排名
  hdac9_rank <- which(markers_sorted$gene == "HDAC9")
  
  # 输出排名
  hdac9_rank
}


{
  load("./subannoComplete.Rdata")
  
  # 提取 meta data
  meta_data <- subsce@meta.data
  
  # 计算 group 中 sdd 和 mdd 在 subcelltype 各类细胞中的占比
  prop_table_1 <- meta_data %>%
    group_by(subcelltype, group) %>%
    summarise(count = n()) %>%
    group_by(subcelltype) %>%
    mutate(prop = count / sum(count))
  
  # 计算 subcelltype 在 group 的 sdd 和 mdd 中的占比
  prop_table_2 <- meta_data %>%
    group_by(group, subcelltype) %>%
    summarise(count = n()) %>%
    group_by(group) %>%
    mutate(prop = count / sum(count))
  
  # 图1：group中sdd和mdd在subcelltype每一类细胞中的占比
  ggplot(prop_table_1, aes(x = subcelltype, y = prop, fill = group)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(title = "Proportion of sdd and mdd in each subcelltype",
         y = "Proportion", x = "Subcelltype") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 图2：subcelltype在group的sdd和mdd中的占比
  ggplot(prop_table_2, aes(x = group, y = prop, fill = subcelltype)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(title = "Proportion of subcelltypes in sdd and mdd groups",
         y = "Proportion", x = "Group") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 图3: SDD和MDD在样本中的umap图
  DimPlot(subsce,reduction = "umap", group.by = "group")
}

{
  #SERPINA1的umap图
  load("./subannoComplete.Rdata")
  
  # 假设你的 Seurat 对象名为 `seurat_object`
  # 确认基因是否存在于数据中
  gene_of_interest <- "SERPINA1"  # 替换为你感兴趣的基因名称
  if (!gene_of_interest %in% rownames(subsce)) {
    stop("基因不存在于Seurat对象中，请检查基因名称。")
  }
  # 查看基因表达值范围
  range(subsce@assays$RNA$data[gene_of_interest, ])
  # 绘制基因表达量的 UMAP 图
  p <- FeaturePlot(
    object = subsce,
    features = gene_of_interest,
    reduction = "umap",
    cols = c("lightgrey" ,"#DE1F1F"),
    #pt.size = 0.5,
    #min.cutoff = 0,
    #max.cutoff = 5
  )
  
  ggsave("./图/SERPINA1的UMAP图.pdf", plot = p, width = 10, height = 10, dpi = 300)
  
}

{
  # 加载必要包
  library(Seurat)
  library(ggplot2)
  library(ggpubr)  # 用于添加显著性标记
  
  # Seurat 对象为 subsce，分组信息和亚群信息存储在 "group" 和 "subcelltype" 列
  load("./subannoComplete.Rdata")
  # 确保基因存在于 Seurat 对象中
  gene <- "SERPINA1"
  if (!(gene %in% rownames(subsce))) {
    stop(paste("基因", gene, "不在 Seurat 对象中，请检查基因名是否正确"))
  }
  
  # 提取表达数据
  expr_data <- FetchData(subsce, vars = c(gene, "group", "subcelltype"))
  
  # 1. 绘制 MDD 和 SDD 中 SERPINA1 表达的箱线图
  # 统计每个样本组的表达
  expr_MDD_SDD <- expr_data[, c(gene, "group")]
  
  # 绘制箱线图
  p1 <- ggplot(expr_MDD_SDD, aes(x = group, y = SERPINA1, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
    #geom_jitter(position = position_jitter(0.2), size = 0.5, alpha = 0.8) +
    scale_fill_manual(values = c("blue", "red")) +  # 自定义颜色
    labs(x = "Sample Type", y = "Expression", title = "SERPINA1 Expression in MDD and SDD") +
    theme_minimal() +
    stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test")  # 添加显著性标记
  
  # 保存图片
  ggsave("./图/SERPINA1_MDD_SDD_Boxplot.pdf", plot = p1, width = 6, height = 4, dpi = 300)
  
  # 2. 绘制 6 个细胞亚群中 SERPINA1 表达的箱线图
  expr_celltypes <- expr_data[, c(gene, "subcelltype")]
  
  # 绘制箱线图
  p2 <- ggplot(expr_celltypes, aes(x = subcelltype, y = SERPINA1, fill = subcelltype)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
    #geom_jitter(position = position_jitter(0.2), size = 0.5, alpha = 0.8) +
    scale_fill_brewer(palette = "Set3") +  # 使用 RColorBrewer 的调色板
    labs(x = "Cell Subtype", y = "Expression", title = "SERPINA1 Expression in 6 Cell Subtypes") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 旋转 x 轴标签
  
  # 保存图片
  ggsave("./图/SERPINA1_Cell_Subtypes_Boxplot.pdf", plot = p2, width = 8, height = 4, dpi = 300)
}