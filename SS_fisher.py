import pandas as pd
import math
import scipy.stats as stats
import sys
import numpy as np

def fisher(set1, set2, bigset1, bigset2):
    # 计算 b, s, p, n
    b_list = find_overlap(set1, set2)
    b = len(b_list)
    all_list = find_union(bigset1, bigset2)
    s = len(set1) - len(b_list)
    p = len(set2) - len(b_list)
    n = len(all_list) - s - p - b
    print(b,s,p,n)
    # 确保 b, s, p 不为 0
    b = b if b != 0 else 1e-10
    s = s if s != 0 else 1e-10
    p = p if p != 0 else 1e-10

    # 构建正确的列联表
    table =[[b,p], [b+s,p+n]]
    #table = [[b, p], [s, n]]

    # 使用 Fisher 精确检验
    fisher_statistic, p_value = stats.fisher_exact(table, alternative='greater')
    if fisher_statistic == 0 or np.isnan(fisher_statistic) or np.isinf(fisher_statistic):
        fisher_statistic = 0.01

    # 根据文献公式计算 SS_Fisher
    # score = math.log((b * n) / (s * p))
    return (fisher_statistic, p_value)

def find_overlap(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    return list(set1 & set2)

def find_union(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    return list(set1 | set2)

#def scale_values(series):
#    """Scales the values of a pandas Series to the range [-1, 1]."""
#    min_val = series.min()
#    max_val = series.max()
#    scaled_series = 2 * ((series - min_val) / (max_val - min_val)) - 1
#    return scaled_series

# 读取输入文件
Lfile = sys.argv[1]
Sfile = sys.argv[2]

# 读入drug gene set L, 找到上调及下调基因
L_df = pd.read_csv(Lfile, header=None, sep='\t')
L_df[1] = pd.to_numeric(L_df[1], errors='coerce').fillna(0).astype(float)

# 缩放药物靶点数据
# L_df[1] = scale_values(L_df[1])

# 定义上调和下调基因的阈值
Lupsig_df = L_df[L_df[1] > 0.5]  # 上调基因
Lup_df = L_df[L_df[1] > 0]  # 上调基因
Ldownsig_df = L_df[L_df[1] < -0.5]  # 下调基因
Ldown_df = L_df[L_df[1] < 0]  # 下调基因

# 读入signture，找到上调及下调基因
S_df = pd.read_csv(Sfile, sep='\t')
S_df['logFC'] = pd.to_numeric(S_df['logFC'], errors='coerce').fillna(0).astype(float)
S_df['pv'] = pd.to_numeric(S_df['pv'], errors='coerce').fillna(0).astype(float)

# 缩放疾病反应数据
# S_df['logFC'] = scale_values(S_df['logFC'])

# 定义上调和下调基因的阈值
Sup_df = S_df[S_df['logFC'] > 0]
Supsig_df = S_df[(S_df['logFC'] > 0.2)]
Sdown_df = S_df[S_df['logFC'] < 0]
Sdownsig_df = S_df[(S_df['logFC'] < -0.2)]

(uu_ratio,uu_pv) = fisher(Supsig_df['Gene_Symbol'], Lupsig_df[0], Sup_df['Gene_Symbol'], Lup_df[0])
(dd_ratio,dd_pv) = fisher(Sdownsig_df['Gene_Symbol'], Ldownsig_df[0], Sdown_df['Gene_Symbol'], Ldown_df[0])
(ud_ratio,ud_pv) = fisher(Supsig_df['Gene_Symbol'], Ldownsig_df[0], Sup_df['Gene_Symbol'], Ldown_df[0])
(du_ratio,du_pv) = fisher(Sdownsig_df['Gene_Symbol'], Lupsig_df[0], Sdown_df['Gene_Symbol'], Lup_df[0])
print(uu_ratio,dd_ratio,ud_ratio,du_ratio)
pro=(math.log(uu_ratio)+math.log(dd_ratio))/2
supp=(math.log(ud_ratio)+math.log(du_ratio))/2 #以e为底

if pro > supp:
    fisher_score = abs(pro)
else:
    fisher_score = -abs(supp)

if fisher_score == np.inf:
    fisher_score = 1
elif fisher_score == -np.inf:
    fisher_score = -1

with open('out_fisher.txt', 'w') as file:
    file.write(str(fisher_score))
