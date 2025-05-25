#计算S（疾病signature）与L（药物靶点）重叠基因的表达水平相关性，两种算法 pearson spearman
#返回 相关系数与p-value
import pandas as pd
import numpy as np
import sys
from scipy.stats import pearsonr
from scipy.stats import spearmanr

def find_overlap(list1, list2):
    # 使用set去除重复元素并找出交集
    set1 = set(list1)
    set2 = set(list2)
    return list(set1 & set2)

Lfile=sys.argv[1]#Lset.txt
Sfile=sys.argv[2]#sig.txt

#读入Disease signture
S_df = pd.read_csv(Sfile, sep='\t')
S_df['logFC'] = pd.to_numeric(S_df['logFC'], errors='coerce').fillna(0).astype(float)#将logfc列变为数值型
#S_df['Gene_Symbol']=S_df['Gene_Symbol'].str.title()#将疾病基因名变为首字母大写
#读入drug gene set L
L_df = pd.read_csv(Lfile, header=None, sep='\t')
L_df[1] = pd.to_numeric(L_df[1], errors='coerce').fillna(0).astype(float)#将药物靶点的权重列变为数值型
#L_df[0]=L_df[0].str.title()#将疾病基因名变为首字母大写

#转换为list
set1 = S_df['Gene_Symbol'].tolist()
set2 = L_df[0].tolist()

ol_list=find_overlap(set1,set2)

#设定gene_symbol为dataframe索引
S_df = S_df.set_index('Gene_Symbol')

#设定第一列为dataframe索引
L_df = L_df.set_index(0)

#筛选在overlap中的S
selected_S_df = S_df.query('index in @ol_list')

#筛选在overlap中的L
selected_L_df = L_df.query('index in @ol_list')

x=list()
y=list()
for index, row in selected_S_df.iterrows():
    if (isinstance(selected_S_df.loc[index]['logFC'], float))&(isinstance(selected_L_df.loc[index][1], float)):
        x.append(selected_S_df.loc[index]['logFC'])
        y.append(selected_L_df.loc[index][1])

pearson_corr, pearson_pv = pearsonr(x, y)
#print(pearson_corr, pearson_pv)

spearman_corr, pearson_pv = spearmanr(x, y)
#print(spearman_corr, pearson_pv)

with open('out_corr.txt', 'w') as file:
    file.write(str(pearson_corr))
    file.write('\n')
    file.write(str(spearman_corr))