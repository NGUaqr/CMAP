import os

# 设置目标文件夹路径
folder_path = r"C://Users//huihui1126//Desktop//药物预测方法测试//cmap//测试数据"  # ← 修改为你自己的路径

# 遍历文件夹中的所有文件
for filename in os.listdir(folder_path):
    if filename.startswith("diseaselist_") and filename.endswith(".txt"):
        old_path = os.path.join(folder_path, filename)
        new_filename = filename.replace("diseaselist_", "druglist_", 1)
        new_path = os.path.join(folder_path, new_filename)
        
        os.rename(old_path, new_path)
        print(f"重命名: {filename} → {new_filename}")
