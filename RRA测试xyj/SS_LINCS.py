#计算疾病的上调基因在药物profile中的GSEA富集情况
import pandas as pd
from gseapy import GSEA
import gseapy as gp
import sys

Lfile=sys.argv[1]#lincs drug profile
Sfile=sys.argv[2]#disease signature
# 测试用 added by zhangbo
# Lfile="C:\\Users\huihui1126\Desktop\cmap_zhangbo\druglist-pos\Acetylsalicylic_acid.txt "#lincs drug profile
# Sfile="C:\\Users\huihui1126\Desktop\cmap_zhangbo\diseaselist-deg2funciton\Granulocyte_DNB_early_before.txt"       #disease signature
# Lfile="C:\\Users\huihui1126\Desktop\cmap-并行化分析\druglist-pos\Acetylsalicylic_acid.txt "#lincs drug profile
# Sfile="C:\\Users\huihui1126\Desktop\cmap-并行化分析\diseaselist-deg2funciton\GSE104187_1d-angionesis.txt"       #disease signature
#读入signture，并找到上调及下调基因
S_df = pd.read_csv(Sfile, sep='\t')
#S_df['Gene_Symbol']=S_df['Gene_Symbol'].str.title()  # 首字母大写 已完成鼠源人源基因转换，注释掉
S_df['logFC'] = pd.to_numeric(S_df['logFC'], errors='coerce').fillna(0).astype(float)
S_df['pv'] = pd.to_numeric(S_df['pv'], errors='coerce').fillna(0).astype(float)
Sup_df = S_df[(S_df['logFC']>0)& (S_df['pv']<0.05)]
Sdown_df = S_df[(S_df['logFC']<-0)& (S_df['pv']<0.05)]
#print(Sup_df["Gene_Symbol"])
#print(Sdown_df["Gene_Symbol"])
with open('gene_sets-defined.gmt', 'w', encoding='utf-8') as file:
    file.write('down'+'\t')
    for item2 in Sdown_df['Gene_Symbol'].tolist():
        file.write('\t'+item2)
    file.write('\n'+'up'+'\t')
    for item1 in Sup_df['Gene_Symbol'].tolist():
        file.write('\t'+item1)
    file.write('\n')
#读入gene set L,并去掉其空值行
L_df = pd.read_csv(Lfile, header=None, sep='\t')
#L_df[0]=L_df[0].str.title() #将药物基因名变为首字母大写 # 首字母大写 已完成鼠源人源基因转换，注释掉
L_df[1] = pd.to_numeric(L_df[1], errors='coerce').fillna(0).astype(float)
L_df_sorted = L_df.sort_values(by=1, ascending=True,ignore_index=False)#根据corr列升序排列
L_df_sorted.dropna(axis=0, how='any', inplace=True)
#print(L_df_sorted)
intersection = list(set(S_df["Gene_Symbol"]) & set(L_df_sorted[0]))
#print(intersection)
pre_res = gp.prerank(rnk=L_df_sorted, gene_sets='gene_sets-defined.gmt',
                     threads=4,min_size=0,max_size=7000,permutation_num=50,
                     outdir='test', format='png', seed=6)
#print(pre_res.res2d.head())




#输出结果
#terms = pre_res.res2d.Term
#由于DNB数量有限,可能只有上调或下调的基因,up或者down可能不存在对应的key,如果没有则赋给对应的ES得分为0
ES_down = pre_res.results.get('down', {}).get('es', 0)
NES_down = pre_res.results.get('down', {}).get('nes', 0)
#print(ES_down, NES_down)  # modified by zhangbo

ES_up = pre_res.results.get('up', {}).get('es', 0)
NES_up = pre_res.results.get('up', {}).get('nes', 0)
#print(ES_up, NES_up)  # modified by zhangbo


# ES_down=pre_res.results['down'].get('es')
# NES_down=pre_res.results['down'].get('nes')
# print(ES_down,NES_down)
# ES_up=pre_res.results['up'].get('es')
# NES_up=pre_res.results['up'].get('nes')  writed by Guo
# 计算综合得分
SSlincs=(ES_up-ES_down)/2
SSlincs_Norm=(NES_up-NES_down)/2
#print(SSlincs,SSlincs_Norm)
#写入文档
with open('out_lincs.txt', 'w') as file:
    file.write(str(SSlincs))
    file.write('\n')
    file.write(str(SSlincs_Norm))