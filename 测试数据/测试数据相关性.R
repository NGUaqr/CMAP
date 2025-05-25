library(dplyr)

# 设置文件夹路径
folder_path <- "C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/测试数据/"

# 获取所有txt文件的文件名
file_list <- list.files(folder_path, pattern = "*.txt", full.names = TRUE)

# 初始化一个空列表来存储数据框
data_list <- list()

# 循环读取每个文件并存储在列表中
for (file in file_list) {
  # 读取文件，假设文件没有列名，所以 use.names = FALSE
  temp_data <- read.table(file, header = FALSE)
  # 提取行名和数值列
  temp_data <- temp_data[, c(1, 2)]
  colnames(temp_data) <- c("RowName", tools::file_path_sans_ext(basename(file)))
  # 添加到列表中
  data_list[[file]] <- temp_data
}

# 合并所有数据框，按行名对齐
merged_data <- Reduce(function(x, y) merge(x, y, by = "RowName", all = TRUE), data_list)

# 将第一列设为行名
rownames(merged_data) <- merged_data$RowName
merged_data <- merged_data[, -1]

# 确保所有数据列都是数值型
merged_data[] <- lapply(merged_data, function(x) as.numeric(as.character(x)))

# 打印合并后的数据框
print(merged_data)
cor_results <- cor(merged_data, use = "pairwise.complete.obs")
# 绘制相关性热图
library(corrplot)
library(RColorBrewer) # 用来配色的
# 挑选相关系数大于阈值的特征对
# 绘制相关热图
corrplot(
  cor_results,
  type = "full",
  # type = c("full", "lower", "upper"),
  #order = "original", #order = c("original", "AOE", "FPC", "hclust", "alphabet"),
  # hclust.method = c("complete", "ward", "ward.D", "ward.D2", "single", "average","mcquitty", "median", "centroid"),
  #tl.col = "black",
  #tl.srt = 45,
  #addCoef.col = "black",
  col = colorRampPalette(c("blue", "red"))(4),
  #number.cex = 0.6,
  #title = basename(input_path[i]) %>% gsub("\\.txt$", "", .),
  #col.lim = c(-1, 1),
  #mar = c(3, 1, 2, 1)
)
dev.off()
# 使用pheatmap绘制热图
pheatmap(cor_results, 
         fontsize = 10, 
         fontsize_number = 8, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_rows = FALSE, 
         cluster_cols = FALSE)