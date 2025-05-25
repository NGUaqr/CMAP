# 常见的小鼠基因列表（这是一个简短的示例，实际使用时请替换为更全面的基因列表）
common_mouse_genes <- c("Akt1", "Bcl2", "Cdkn1a", "Foxp3", "Gata3", 
                        "Il10", "Il6", "Mapk1", "Pten", "Tnf", 
                        "Ccl2", "Ccl5", "Cxcl10", "Cxcr4", "Egr1", 
                        "Fos", "Jun", "Mtor", "Nfkb1", "Stat3")

# 设置随机种子以保证结果可重复
set.seed(123)

# 随机选择20个基因
selected_genes <- sample(common_mouse_genes, 20)

# 为这些基因生成随机的log2FC值
log2fc_values <- runif(20, min = -2, max = 2)  # 随机生成[-2, 2]之间的值

# 按照log2FC值从大到小排序
sorted_indices <- order(-log2fc_values)
selected_genes <- selected_genes[sorted_indices]
log2fc_values <- log2fc_values[sorted_indices]

# 保证前10个log2FC为正，后10个为负
log2fc_values[1:10] <- abs(log2fc_values[1:10])
log2fc_values[11:20] <- -abs(log2fc_values[11:20])

# 构建疾病数据框
disease_df <- data.frame(Gene = selected_genes, log2FC = log2fc_values)

# 保存疾病文件
write.table(disease_df, file = "disease_file.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# 生成20个药物文件
for (i in 0:19) {
  drug_log2fc <- log2fc_values
  if (i > 0) {
    # 将最小的i个log2FC符号修改为正（从完全一致到逐步反转）
    drug_log2fc[(20-i+1):20] <- -drug_log2fc[(20-i+1):20]
  }
  drug_df <- data.frame(Gene = selected_genes, log2FC = drug_log2fc)
  drug_file_name <- paste0("drug_file_", sprintf("%02d", i+1), ".txt")
  write.table(drug_df, file = drug_file_name, row.names = FALSE, col.names = F,sep = "\t", quote = FALSE)
}

cat("Disease and drug files generated successfully!\n")
