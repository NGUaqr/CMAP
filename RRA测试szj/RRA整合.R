library(RobustRankAggreg)
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(stringr)
library(pheatmap)
library(coin)
library(readr)
install.packages("extrafont")
library(extrafont)
### 正效应排序-得分越大正效应越强,加个负号,排名越靠前
input_path <- "C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/RRA测试szj/output/正效应/"
output_path <- "C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/RRA测试szj/output/"
filename <- list.files(input_path)
input_path <- paste(input_path,filename,sep = '')
#读入result.txt到dataframe
i=1
for (i in 1:length(input_path)) {
  print(input_path[i])
  new_df <- df <-read.table(file = input_path[i],row.names = 1,sep="\t",header=TRUE, quote = "", fileEncoding = "GBK")
  new_df<- df <- df[,-c(1,2)]### 去掉fisher和cmap
  for (col_name in names(df)) {  
    if(grepl("SS..", col_name)){
      #加个负号,排名靠前的是正效应越强
      new_df[col_name] <- rank(-df[col_name])/length(order(df[[col_name]]))	
    }
    
  }
  # 将数据帧转换为矩阵
  mat <- as.matrix(new_df)
  #print(mat)
  #计算稳健整合分数
  # write.table(new_df,file = paste(output_path, "rank", filename[i], sep = ''), row.names = TRUE, col.names = TRUE, sep = "\t", quote = TRUE) 
  ranksum<- aggregateRanks(rmat = mat, method = "RRA")
  write.table(ranksum,file = paste(output_path,filename[i],sep = ''), row.names = TRUE, col.names = TRUE, sep = "\t", quote = TRUE) 
} 
#### 绘制排序热图和RAA整合热图
pdf("C:/Users/huihui1126/Desktop/5vs11指标原始数据.pdf", width = 8, height = 66, family="GB1")
pheatmap(df, scale = "column",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         cluster_rows = T, 
         cluster_cols = FALSE)

dev.off()
pdf("C:/Users/huihui1126/Desktop/5vs11指标原始数据加负号排序.pdf", width = 8, height = 66, family="GB1")
pheatmap(mat, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         cluster_rows = T, 
         cluster_cols = FALSE)
write.csv(mat,"2vs11指标原始数据加负号排序.csv")
dev.off()
pdf("C:/Users/huihui1126/Desktop/5vs11指标RRA整合排序.pdf", width = 8, height = 66, family="GB1")
pheatmap(ranksum[,-1,drop=FALSE], 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         cluster_rows = F, 
         cluster_cols = FALSE)
dev.off()
### 负效应排序-得分越小负效应越强,排名越靠前
input_path <- "C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/RRA测试szj/output/负效应/"
output_path <- "C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/RRA测试szj/output/"
filename <- list.files(input_path)
input_path <- paste(input_path,filename,sep = '')
#读入result.txt到dataframe
i=1
for (i in 1:length(input_path)) {
  print(input_path[i])
  new_df<-df<-read.table(file = input_path[i],row.names = 1,sep="\t",header=TRUE, quote = "", fileEncoding = "GBK")
  new_df<- df <- df[,-c(1,2)]### 去掉fisher和cmap
  for (col_name in names(df)) { 
    if(grepl("SS..", col_name)){
      #print(df[[col_name]])
      #排名靠前的是负效应越强
      new_df[col_name] <- rank(df[col_name])/length(order(df[[col_name]]))
      #print(df[col_name])
    }
    
  }
  # 将数据帧转换为矩阵
  mat <- as.matrix(new_df)
  #print(mat)
  #计算稳健整合分数
  # write.table(new_df,file = paste(output_path, "rank", filename[i], sep = ''), row.names = TRUE, col.names = TRUE, sep = "\t", quote = TRUE) 
  ranksum<- aggregateRanks(rmat = mat, method = "RRA")
  write.table(ranksum,file = paste(output_path,filename[i],sep = ''), row.names = TRUE, col.names = TRUE, sep = "\t", quote = TRUE) 
} 
#### 绘制排序热图和RAA整合热图
pheatmap(df, scale = "column",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         cluster_rows = T, 
         cluster_cols = FALSE)
pheatmap(mat, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         cluster_rows = T, 
         cluster_cols = FALSE)
pheatmap(ranksum[,-1,drop=FALSE], 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         cluster_rows = F, 
         cluster_cols = FALSE)

#### 整合before与after ####
dev.off()
# 读取药物通路数据和时间注释
#moa <- read.table("C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/RRA测试/moa.txt", header = TRUE, row.names = 1, sep = "\t")
moa <- read.table("C:/Users/huihui1126/Desktop/药物预测方法测试/moa.txt", header = TRUE, row.names = 1, sep = "\t")
moa <- read.table("C:/Users/huihui1126/Desktop/药物预测方法测试/moa103.txt", header = TRUE, row.names = 1, sep = "\t")
annotation_row <- as.data.frame(moa)
output_path <- "C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/RRA测试/output/"
filename  <- list.files(output_path,pattern = "\\.txt$")
filename <- paste(output_path,filename,sep = '')
#### before和after单边测试
i=1
for (i in 1:length(filename)) {
  annotation_row <- as.data.frame(moa)
  output_path <- "C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/RRA测试/output/"
  print(filename[i])
  df <- read.table(file = filename[i],row.names = 1,sep="\t",header=TRUE)
  df <- subset(df,Score<mean(Score)) 
  # 将 Name 列设置为行名并删除 .txt 后缀
  rownames(df) <- gsub("\\.txt$", "", df$Name)
  df <- df[, -1, drop = FALSE]  # 保持数据框结构
  rownames(df) <- gsub("-", "_", rownames(df))
  df$Score <- rank(df$Score, ties.method = "average")
  Pathway_color <- c(
    "MAPK" = "red", 
    "PI3K/Akt/Nrf2/HO-1" = "blue", 
    "RhoA/ROCK" = "green", 
    "Sonic Hedgehog" = "purple", 
    "Hippo/YAP" = "orange", 
    "NF-KB" = "brown", 
    "Wnt/b-catenin" = "pink", 
    "TLR4" = "yellow", 
    "NLRP3/IL-1b" = "cyan", 
    "Sildenafil" = "#5F80B4", 
    "TGF-b/SMADs" = "black", 
    "NRF2/HO-1" = "violet", 
    "PI3K/Akt/mTOR" = "magenta", 
    "TGFb1/TAK1,TGF-b2,TGF-b3,RhoA/ROCK" = "navy",
    "NA" = "white")
  Time_color <- c("early" ="red", "late" = "blue", "long_term" = "green","NA" = "white")
  ann_colors <- list(#Pathway=Pathway_color, 
                     Time= Time_color) #颜色设置
  # 确保行注释的顺序与矩阵行的顺序一致
  annotation_row <- annotation_row[rownames(df), ]
  all(rownames(df) == rownames(annotation_row)) # 应返回TRUE
  # 绘制热图
  pheatmap(as.matrix(df),
           scale = "none",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           color=colorRampPalette(c("firebrick3","white","navy"))(10),
           annotation_row = annotation_row,
           annotation_colors = ann_colors,
           main = basename(filename[i]),
           fontsize = 10)
  #提取early类药物和其他药物的排名进行比较：
  early_ranks <- df$Score[annotation_row$Time %in% c('early',"long_term")]
  other_ranks <- df$Score[annotation_row$Time %in% c('late',"long_term")]
  early_mean_rank <- mean(early_ranks)
  other_mean_rank <- mean(other_ranks)
  relative_percentage <- (early_mean_rank / other_mean_rank) * 100
  print(paste("Early类药物的平均排名得分:", early_mean_rank))
  print(paste("late类药物的平均排名得分:", other_mean_rank))
  print(paste("Early-late相对百分比:", relative_percentage))
  #进行Mann-Whitney U检验：
  df$Category <- annotation_row$Time
  w_df <- df[df$Category != 'long_term', ]
  w_df$Category <- as.factor(w_df$Category)
  test_result <- wilcox_test(Score~Category, data = w_df, distribution = "exact")
  p_value <- pvalue(test_result)
  print(paste("Mann-Whitney U检验 p值:", p_value))
  #提取实验证据类药物和无实验证据药物的排名进行比较：
  clinical_ranks <-  df[!is.na(annotation_row$Pathway),] %>% .[["Score"]]
  NA_ranks <-  df[is.na(annotation_row$Pathway),] %>% .[["Score"]]
  clinical_mean_ranks <- mean(clinical_ranks)
  NA_mean_ranks <- mean(NA_ranks)
  relative_percentage <- (clinical_mean_ranks / NA_mean_ranks) * 100
  print(paste("实验证据类药物的平均排名得分:", clinical_mean_ranks))
  print(paste("无证据类药物的平均排名得分:", NA_mean_ranks))
  print(paste("实验-无证据相对百分比:", relative_percentage))
  #进行Mann-Whitney U检验：
  df$Category <- annotation_row$Pathway
  df[!is.na(df$Category),][["Category"]] <- 'clinical'
  df[is.na(df$Category),][["Category"]] <- '无'
  df$Category <- as.factor(df$Category)
  test_result <- wilcox_test(Score~Category, data = df, distribution = "exact")
  p_value <- pvalue(test_result)
  print(paste("Mann-Whitney U检验 p值:", p_value))
  
}

dev.off()


#### 合并取min排序作为综合打分 ####
i = 1
for (i in seq(1, length(filename), 2)) {
  annotation_row <- as.data.frame(moa)
  output_path <- "C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/RRA测试/output/"
  # 读取第 i 个文件
  file1 <- filename[i]
  data1 <- read.table(file1, header = TRUE, sep = "\t")
  data1 <- subset(data1,Score<mean(Score))
  # 读取第 i+1 个文件
  file2 <- filename[i+1]
  data2 <- read.table(file2, header = TRUE, sep = "\t")
  data2 <- subset(data2,Score<mean(Score))
  #求 Name 列的并集
  common_names <- unique(union(data1$Name, data2$Name))
  # # 根据交集创建一个新的数据框
  merged_data <- data.frame(
    Name = common_names,
    score_1 = data1[match(common_names, data1$Name), "Score"],
    score_2 = data2[match(common_names, data2$Name), "Score"]
  )
  # 按照 score_1 进行由小到大排序,并计算rank
  merged_data <- merged_data[order(merged_data$score_1), ]
  merged_data$rank_score_1 <- rank(merged_data$score_1, ties.method = "average")
  
  # 按照 score_2 进行由小到大排序,并计算rank
  merged_data <- merged_data[order(merged_data$score_2), ]
  merged_data$rank_score_2 <- rank(merged_data$score_2, ties.method = "average")
  # 计算 rank_min
  merged_data$rank_min <- pmin(merged_data$rank_score_1, merged_data$rank_score_2)
  # 按照 rank_avg 进行由小到大排序
  merged_data <- merged_data[order(merged_data$rank_min), ]
  df <- merged_data[,c(1,6)]
  rownames(df) <- gsub("\\.txt$", "", df$Name)
  df <- df[, -1, drop = FALSE]  # 保持数据框结构
  rownames(df) <- gsub("-", "_", rownames(df))
  Pathway_color <- c(
    "MAPK" = "red", 
    "PI3K/Akt/Nrf2/HO-1" = "blue", 
    "RhoA/ROCK" = "green", 
    "Sonic Hedgehog" = "purple", 
    "Hippo/YAP" = "orange", 
    "NF-KB" = "brown", 
    "Wnt/b-catenin" = "pink", 
    "TLR4" = "yellow", 
    "NLRP3/IL-1b" = "cyan", 
    "Sildenafil" = "#5F80B4", 
    "TGF-b/SMADs" = "black", 
    "NRF2/HO-1" = "violet", 
    "PI3K/Akt/mTOR" = "magenta", 
    "TGFb1/TAK1,TGF-b2,TGF-b3,RhoA/ROCK" = "navy",
    "NA" = "white")
  Time_color <- c("early" ="red", "late" = "blue", "long_term" = "green","NA" = "white")
  ann_colors <- list(Pathway=Pathway_color, Time= Time_color) #颜色设置
  # 确保行注释的顺序与矩阵行的顺序一致
  annotation_row <- annotation_row[rownames(df), ]
  all(rownames(df) == rownames(annotation_row)) # 应返回TRUE
  # 绘制热图
  pheatmap(as.matrix(df),
           scale = "none",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           color=colorRampPalette(c("firebrick3","white","navy"))(100),
           annotation_row = annotation_row,
           annotation_colors = ann_colors,
           main = basename(paste0(basename(file1), "_", basename(file2))),
           fontsize = 10)
  #提取early类药物和其他药物的排名进行比较：
  print(paste0(basename(file1), "_", basename(file2)))
  
  early_ranks <- df$rank_min[annotation_row$Time %in% c('early')]
  other_ranks <- df$rank_min[annotation_row$Time %in% c('late')]
  early_mean_rank <- mean(early_ranks)
  other_mean_rank <- mean(other_ranks)
  relative_percentage <- (early_mean_rank / other_mean_rank) * 100
  print(paste("Early类药物的平均排名得分:", early_mean_rank))
  print(paste("late类药物的平均排名得分:", other_mean_rank))
  print(paste("Early-late相对百分比:", relative_percentage))
  #进行Mann-Whitney U检验：
  df$Category <- annotation_row$Time
  w_df <- df
  w_df <- df[df$Category != 'long_term', ]
  w_df$Category <- as.factor(w_df$Category)
  test_result <- wilcox_test(rank_min~Category, data = w_df, distribution = "exact")
  p_value <- pvalue(test_result)
  print(paste("Mann-Whitney U检验 p值:", p_value))
  #提取实验证据类药物和无实验证据药物的排名进行比较：
  clinical_ranks <-  df[!is.na(annotation_row$Pathway),] %>% .[["rank_min"]]
  NA_ranks <-  df[is.na(annotation_row$Pathway),] %>% .[["rank_min"]]
  clinical_mean_ranks <- mean(clinical_ranks)
  NA_mean_ranks <- mean(NA_ranks)
  relative_percentage <- (clinical_mean_ranks / NA_mean_ranks) * 100
  print(paste("实验证据类药物的平均排名得分:", clinical_mean_ranks))
  print(paste("无证据类药物的平均排名得分:", NA_mean_ranks))
  print(paste("实验-无证据相对百分比:", relative_percentage))
  #进行Mann-Whitney U检验：
  df$Category <- annotation_row$Pathway
  df[!is.na(df$Category),][["Category"]] <- 'clinical'
  df[is.na(df$Category),][["Category"]] <- '无'
  df$Category <- as.factor(df$Category)
  test_result <- wilcox_test(rank_min~Category, data = df, distribution = "exact")
  p_value <- pvalue(test_result)
  print(paste("Mann-Whitney U检验 p值:", p_value))
  cat("\n")
}

dev.off()

#### 496药物预测 ####
moa <- read_excel("C:/Users/huihui1126/Desktop/药物预测方法测试/MOA496.xlsx")
moa <- moa[c("中文名","Name")]
output_path <- "C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/RRA测试/output/"
filename  <- list.files(output_path,pattern = "\\.txt$")
filename <- paste(output_path,filename,sep = '')
#### before和after单边测试
i=1
for (i in 1:length(filename)) {
  annotation_row <- as.data.frame(moa)
  output_path <- "C:/Users/huihui1126/Desktop/药物预测方法测试/cmap/RRA测试/output/"
  print(filename[i])
  df <- read.table(file = filename[i],row.names = 1,sep="\t",header=TRUE)
  # 将 Name 列设置为行名并删除 .txt 后缀
  rownames(df) <- gsub("\\.txt$", "", df$Name)
  rownames(df) <- gsub("drug_", "", rownames(df))
  df <- df[, -1, drop = FALSE]  # 保持数据框结构
  rownames(df) <- gsub("-", "_", rownames(df))
  # 1. 将行名转换为普通列
  df$compound_id <- rownames(df)
  annotation_row$compound_id <- rownames(annotation_row)
  # 2. 根据 'compound_id' 进行合并
  merged_df <- merge(df, annotation_row, by = "compound_id", all = TRUE)
  merged_df <- na.omit(merged_df)
  rownames(merged_df) <- merged_df$中文名
  merged_df$Score <- rank(merged_df$Score, ties.method = "average")
  merged_df <- subset(merged_df,Score<=50) 
  merged_df <- merged_df %>% arrange(Score)
  # 绘制热图
  pheatmap(as.matrix(merged_df[,"Score", drop = FALSE]),
           scale = "none",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           color=colorRampPalette(c("firebrick3","white","navy"))(10),
           main = basename(filename[i]),
           fontsize = 10)
}

dev.off()