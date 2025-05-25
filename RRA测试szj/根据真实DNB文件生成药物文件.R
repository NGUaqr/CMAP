setwd("C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/RRA测试/druglist-pos/")
# 读取疾病谱文件
disease_data <- read.table("C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/RRA测试/diseaselist-deg2function/Macrophages_DNB_early_before.txt", header = TRUE, sep = "\t")
disease_data <- read.table("C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/RRA测试/diseaselist-deg2function/Macrophages_DNB_early_after.txt", header = TRUE, sep = "\t")
# 提取靶点名和logFC（活性）
targets <- disease_data$Gene_Symbol
activity <- disease_data$logFC

# 活性的数量
n <- length(activity)

# 生成药物文件
for (i in 1:10) {
  drug_activity <- activity
  
  # 引入反向值
  if (i > 1) {
    indices_to_reverse <- sample(1:n, size = floor(n * (i-1)/10))  # 随机选择一些活性反向
    drug_activity[indices_to_reverse] <- -drug_activity[indices_to_reverse]
  }
  
  # 将活性值放缩到 -1 到 1 之间
  min_val <- min(drug_activity)
  max_val <- max(drug_activity)
  drug_activity <- 2 * ((drug_activity - min_val) / (max_val - min_val)) - 1
  
  # 创建数据框并排序
  drug_data <- data.frame(Target = targets, Activity = drug_activity)
  drug_data <- drug_data[order(-drug_data$Activity), ]  # 按照活性值从大到小排序
  
  # 保存药物文件
  file_name <- paste0("drug_", i+10, ".txt")
  write.table(drug_data, file_name, sep = "\t", col.names = F,row.names = FALSE, quote = FALSE)
}