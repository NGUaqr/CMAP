import pandas as pd
import numpy as np
import statistics
import sys
import re

Lfile = sys.argv[1]
Sfile = sys.argv[2]

# 读入Disease signature
S_df = pd.read_csv(Sfile, sep='\t')
S_df['logFC'] = pd.to_numeric(S_df['logFC'], errors='coerce').fillna(0).astype(float)
S_df = S_df.drop_duplicates(subset='Gene_Symbol')  # 去重

# 读入drug gene set L
L_df = pd.read_csv(Lfile, header=None, sep='\t')
L_df = L_df.drop_duplicates(subset=0)

def zhangscore(S_df, L_df):
    # 转换为list
    set1 = S_df['Gene_Symbol'].tolist()
    set2 = L_df[0].tolist()
    ol_list = find_overlap(set1, set2)

    # 设定gene_symbol为dataframe索引
    S_df = S_df.set_index('Gene_Symbol')
    L_df = L_df.set_index(0)

    # 选择交集基因
    selected_S_df = S_df.loc[ol_list]

    # 按照'logFC'列进行排序，但保持原始索引顺序
    selected_S_df_sorted = selected_S_df.sort_values(by='logFC', ignore_index=False)
    L_df_sorted = L_df.sort_values(by=1, ignore_index=False)

    R_gi = []
    s_gi = []

    # 遍历selected_S_df_sorted，并查找其元素在L_df_sorted的rank
    for index_S in selected_S_df_sorted.index:
        if re.match(r"^[.]$", index_S) is None:
            if index_S in L_df_sorted.index:
                index_position = L_df_sorted.index.get_loc(index_S)
                R_gi.append(index_position + 1)  # Rank starts from 1
                if selected_S_df_sorted.loc[index_S]['logFC'] * L_df_sorted.loc[index_S][1] > 0:
                    s_gi.append(1)
                else:
                    s_gi.append(-1)

    if not R_gi or not s_gi:
        return 0

# 计算公式的分子和分母
    numerator = sum([R_gi[i] * s_gi[i] for i in range(len(R_gi))])
    denominator = sum([(len(R_gi) - i) for i in range(len(R_gi))])
    score = numerator / denominator
    return score

def find_overlap(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    overlap = set1 & set2
    return list(overlap)

# 计算原始得分
S_L_score = zhangscore(S_df, L_df)

#建立随机的drug gene set L_df random 10次
rand=[]
for i in range(10):
    random_L_df=L_df
    random_L_df[0]=np.random.permutation(random_L_df[0])
    rand.append(zhangscore(S_df,random_L_df))
mean_rand=sum(rand) / len(rand)
norm_S_L_score=(S_L_score-mean_rand)/statistics.stdev(rand)
#print(S_L_score,norm_S_L_score)

with open('out_zhang.txt', 'w') as file:
    file.write(f"{S_L_score}\n")
    file.write(f"{norm_S_L_score}\n")
