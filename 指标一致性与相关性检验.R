#### 指标一致性检验 ####
# Load necessary libraries
library(ggplot2)
library(reshape2)
# Define the directory containing the txt files
input_path <- "C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/output/"
files <- list.files(input_path, pattern = "\\.txt$", full.names = TRUE)
data_list <- list()
for (file in files) {
  data <- read.table(file, header = TRUE, sep = "\t", row.names = 1)
  file_name <- basename(file)
  file_name <- sub("\\.txt$", "", file_name)
  data$Filename <- file_name
  data_list[[file_name]] <- data
}

# Combine all data frames into a single data frame
combined_data <- do.call(rbind, data_list)
rownames(combined_data) <- combined_data$Filename
combined_data$Filename <- NULL

# Convert rownames to a column
combined_data <- cbind(Filename = rownames(combined_data), combined_data)
#scale_combined_data <- scale(combined_data[,-1])
# Melt the data frame for easier plotting
#melted_data <- melt(scale_combined_data, id.vars = row.names(combined_data))
melted_data <- melt(combined_data, id.vars = "Filename")
melted_data$Filename <- sub(".*druglist_(\\d+)\\.txt", "\\1", melted_data$Filename)
# Convert Filename to a factor with levels in the correct order
melted_data$Filename <- factor(melted_data$Filename, levels = unique(melted_data$Filename))

# Plot the data
ggplot(melted_data, aes(x = Filename, y = value, group = variable, color = variable)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ variable, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Line Plot for Each Column", x = "Filename", y = "Value") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# Save the plots
ggsave("combined_line_plots.png")
# 生成趋势图
# 替换为你实际的数据列名
melted_data$Index <- as.numeric(sub(".*druglist_(\\d+)\\.txt", "\\1", melted_data$Filename))

# 检查是否提取成功
if (all(is.na(melted_data$Index))) stop("Index column contains only NA. Check Filename format.")
ggplot(melted_data, aes(x = Index, y = value, color = variable, group = variable)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = seq(min(melted_data$Index, na.rm = TRUE),
                                  max(melted_data$Index, na.rm = TRUE),
                                  by = 5)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    title = "Trend of Connectivity Metrics Across Drug Profiles",
    x = "Drug Profile Index",
    y = "Metric Score"
  )


#### 指标相关性检验 ####
rm(list=ls())
## 首先确认下各指标的相关性 以选择相关性最高的指标组合 ####
library(RobustRankAggreg)
library(readxl)
library(gridExtra)
library(psych)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tibble)
folder_path <- "C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/output/"
files <-list.files(path = folder_path,pattern = "\\.txt$",full.names = T)
dataframes <- list()
#i=1
dev.off()
for (file in files){
  print(file)
  df <- read.table(file,sep = "\t",header = T)
  df$filename <- basename(file)
  dataframes[[length(dataframes)+1]] <- df
  }
  combined_df <- bind_rows(dataframes)
  combined_df <- column_to_rownames(combined_df,var = "filename")
  combined_df <- combined_df[,-1]
  print(combined_df)
  corr <- cor(combined_df)
  # 绘制相关性热图
  library(corrplot)
  library(RColorBrewer) # 用来配色的
  # 挑选相关系数大于阈值的特征对
  # 绘制相关热图
  corrplot(
    corr,
    type = "full",
    # type = c("full", "lower", "upper"),
    #order = "original", #order = c("original", "AOE", "FPC", "hclust", "alphabet"),
    # hclust.method = c("complete", "ward", "ward.D", "ward.D2", "single", "average","mcquitty", "median", "centroid"),
    tl.col = "black",
    #tl.srt = 45,
    addCoef.col = "black",
    col = colorRampPalette(c("blue", "red"))(4),
    number.cex = 0.6,
    title = basename(input_path[i]) %>% gsub("\\.txt$", "", .),
    col.lim = c(-1, 1),
    mar = c(3, 1, 2, 1)
  )
  corr <- cor(combined_df, use = "pairwise.complete.obs")
  # 使用pheatmap绘制热图，不进行聚类
  pheatmap(corr, 
           display_numbers = TRUE, 
           fontsize = 10, 
           fontsize_number = 8, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           cluster_rows = FALSE, 
           cluster_cols = FALSE)

## R代码美化##
library(formatR)
# tidy_source("file path")
tidy_source("C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/指标一致性与相关性检验.R")
  


#### 指标一致性检验 ####
library(ggplot2)
library(reshape2)

# 定义数据路径
input_path <- "C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/output/"
files <- list.files(input_path, pattern = "\\.txt$", full.names = TRUE)

# 读取所有文件并加标签
data_list <- list()
for (file in files) {
  data <- read.table(file, header = TRUE, sep = "\t", row.names = 1)
  file_name <- tools::file_path_sans_ext(basename(file))  # 去掉扩展名
  data$Filename <- file_name
  data_list[[file_name]] <- data
}

# 合并数据
combined_data <- do.call(rbind, data_list)
rownames(combined_data) <- combined_data$Filename
combined_data$Filename <- NULL
combined_data <- cbind(Filename = rownames(combined_data), combined_data)

# 提取编号作为分组因子
combined_data$Group <- sub(".*_(\\d+)$", "\\1", combined_data$Filename)

# 重整数据格式
melted_data <- melt(combined_data, id.vars = c("Filename", "Group"))

# 因子化分组变量，保证顺序一致
melted_data$Group <- factor(melted_data$Group, levels = unique(melted_data$Group))

# 绘图
p <- ggplot(melted_data, aes(x = Group, y = value, group = Filename, color = Filename)) +
  geom_line(alpha = 0.6, size = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none",  # 可根据需要保留
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "指标一致性折线图",
    x = "文件编号",
    y = "指标值"
  )

print(p)

