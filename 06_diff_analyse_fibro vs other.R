library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
#加载数据

load("./subannoComplete.Rdata")
# 查看 metadata，确保包含分组信息
head(subsce@meta.data)

# 确保 subcelltype 列中有 fibro 和 其他 组信息
table(subsce@meta.data$subcelltype)


# 找出实验组（SDD）和对照组（MDD）之间的差异基因
diff_genes <- FindMarkers(subsce,ident.1 = "fibro", group.by = "subcelltype",min.pct = 0.1) %>%
  mutate(gene = rownames(.))
# 查看前几个差异基因
head(diff_genes)

# # 从Seurat对象提取基因的平均表达量
# # 假设 seurat_obj 是你的Seurat对象
# avg_exp <- rowMeans(GetAssayData(subsce, layer = "count"))
# 
# # 将平均表达量添加到差异分析结果中
# diff_genes$avg_expression <- avg_exp[rownames(diff_genes)]
# 
# # 筛选表达量大于10的基因
# diff_genes <- diff_genes[diff_genes$avg_expression > 0.1, ]







library(ggplot2)


diff_genes_filter <- diff_genes %>%
  filter(pct.1 > 0.1& p_val_adj < 0.05) %>%
  filter(abs(avg_log2FC) > 0.25)

diff_genes_filter_test <- diff_genes %>%
  filter(pct.1 > 0.1& p_val_adj < 0.05)
diff_genes_filter_test$Delta_Percentage <- (diff_genes_filter_test$pct.1 - diff_genes_filter_test$pct.2)

#top_left <- diff_genes_filter_test %>%
#  arrange(avg_log2FC, Delta_Percentage) %>%
#  slice_head(n = 500)

#####
# Min-Max 归一化函数
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

# 对 logFC 和 Delta_Percentage 进行归一化
diff_genes_filter_test <- diff_genes_filter_test %>%
  mutate(
    normalized_logFC = normalize(avg_log2FC),
    normalized_Delta_Percentage = normalize(Delta_Percentage)
  )

# 加权综合考虑 80% 的 normalized_logFC 和 20% 的 normalized_Delta_Percentage
weight_logFC <- 0.5
weight_Delta_Percentage <- 0.5

diff_genes_filter_test <- diff_genes_filter_test %>%
  mutate(weighted_score = weight_logFC * normalized_logFC + weight_Delta_Percentage * normalized_Delta_Percentage)

# 按照加权分数排序并获取前 100 个
top_left <- diff_genes_filter_test %>%
  arrange(weighted_score) %>%  # 排序
  slice_head(n = 100)  # 取前 100 个
####
#存储topleft信息
save(top_left,file = "./fibrovsother差异基因.Rdata")
# 设定需要标红的基因列表
highlight_genes <- rownames(top_left)  # 替换为实际的基因名称
diff_genes_filter_test$gene_name <- rownames(diff_genes_filter_test)

# 添加一个列以确定基因是否在高亮列表中
diff_genes_filter_test$highlight <- ifelse(diff_genes_filter_test$gene_name %in% highlight_genes, "highlight", "normal")

# 绘制散点图，标红高亮基因
ggplot(diff_genes_filter_test, aes(x = avg_log2FC, y = Delta_Percentage, color = highlight)) +
  geom_point() +
  scale_color_manual(values = c("highlight" = "red", "normal" = "blue")) +  # 高亮基因标红，其他标蓝
  theme_minimal() +
  labs(title = "Scatter Plot with Highlighted Genes",
       x = "avg_log2FC",
       y = "Delta_Percentage",
       color = "Gene Status")
# 保存前 100 个基因到 CSV 文件
write.csv(highlight_genes, file = "差异基因fibro vs other.csv", row.names = FALSE)


#end of current process
#############################################################

library(ggrepel)
colnames(diff_genes_filter)
P.Value_t = 1e-38
# 添加log fold change和p值阈值的列，用于绘图
diff_genes_filter$logFC <- diff_genes_filter$avg_log2FC
diff_genes_filter$Significance <- ifelse(diff_genes_filter$p_val_adj < 0.05 & abs(diff_genes_filter$logFC) > 0.5, "Significant", "Not Significant")
table(diff_genes_filter$Significance)
# 绘制火山图
ggplot(diff_genes_filter, aes(x = logFC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano plot", x = "Log Fold Change", y = "-Log10 Adjusted P-Value") +
  geom_hline(yintercept = -log10(P.Value_t)) + 
  theme(legend.position = "right")



# 筛选显著差异表达的基因，p值 < 0.05 且 logFC > 0.25
sig_genes <- rownames(diff_genes_filter[diff_genes_filter$p_val_adj < 0.01 & abs(diff_genes_filter$logFC) > 0.25, ])

####测试percentage difference
diff_genes_filter$Delta_Percentage <- (diff_genes_filter$pct.1 - diff_genes_filter$pct.2)
plot(diff_genes_filter$logFC,diff_genes_filter$Delta_Percentage)


# 根据 p 值从小到大排序
#sorted_diff_genes <- diff_genes_filter[order(diff_genes_filter$p_val_adj), ]

# 提取前 100 个基因
#top_100_genes <- rownames(head(sorted_diff_genes, 300))
#提取下调的前N个基因
topn_downregulated <- diff_genes_filter %>%
  filter(avg_log2FC < 0) %>%  # 筛选下调基因
  arrange(avg_log2FC) %>%  # 按log2FC从小到大排序
  slice_head(n = 300)  # 取前10个基因

topn_downregulated <- rownames(topn_downregulated)
write.csv(topn_downregulated, file = "provsfibro_top_100_differential_genes.csv", row.names = FALSE)








#表达量上调和下调的
# 筛选出上调的前10个基因
top_upregulated <- diff_genes_filter %>%
  filter(avg_log2FC > 0) %>%  # 筛选上调基因
  arrange(desc(avg_log2FC)) %>%  # 按log2FC从大到小排序
  slice_head(n = 10)  # 取前10个基因

# 筛选出下调的前10个基因
top_downregulated <- diff_genes_filter %>%
  filter(avg_log2FC < 0) %>%  # 筛选下调基因
  arrange(avg_log2FC) %>%  # 按log2FC从小到大排序
  slice_head(n = 10)  # 取前10个基因
# 为上调和下调基因名称添加标签
top_upregulated$gene_name <- paste("[Upregulated]", rownames(top_upregulated))
top_downregulated$gene_name <- paste("[Downregulated]", rownames(top_downregulated))

# 合并上调和下调的基因列表到数据框
top_genes_df <- data.frame(
  gene_name = c(top_upregulated$gene_name, top_downregulated$gene_name)
)
# 合并上调和下调基因
top_genes <- c(rownames(top_upregulated), rownames(top_downregulated))
# 将基因保存为CSV文件
write.csv(top_genes_df, file = "FIBROvsPRO_top_20_genes.csv", row.names = FALSE)

# 使用 Seurat 的 DoHeatmap 函数绘制热图
#画热图的时候提示某些基因在缩放的数据里找不到，而在前面缩放默认缩放的是找到的前2000个高变基因

Idents(subsce) <- "subcelltype"
 DoHeatmap(subsce, features = top_genes,group.by = "subcelltype") +
   scale_fill_gradientn(colors = c("blue", "white", "red")) +
  labs(title = "Differentially Expressed Genes Heatmap")

# 查看前 100 个基因
print(top_100_genes)

# 使用 Seurat 绘制前 100 个差异基因的热图
# DoHeatmap(subsce, features = top_100_genes) +
#   scale_fill_gradientn(colors = c("blue", "white", "red")) +
#   labs(title = "Top 100 Differentially Expressed Genes Heatmap")

# 保存前 100 个基因到 CSV 文件
write.csv(top_100_genes, file = "provsfibro_top_100_differential_genes.csv", row.names = FALSE)
