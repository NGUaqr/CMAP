
import os

# 设置目标文件夹路径
folder_path = r"C:\Users\huihui1126\Desktop\药物预测方法测试\cmap\druglist-pos"  # ← 替换为你的实际路径

# 遍历所有文件
for filename in os.listdir(folder_path):
    if filename.startswith("druglist_") and filename.endswith(".txt"):
        file_path = os.path.join(folder_path, filename)

        new_lines = []
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split()  # 可改为 split('\t') 如果是制表符分隔
                if len(parts) >= 3:
                    # 删除第三列
                    parts.pop(2)
                new_lines.append('\t'.join(parts) + '\n')  # 使用 \t 拼接，可改为空格

        with open(file_path, 'w', encoding='utf-8') as f:
            f.writelines(new_lines)

        print(f"已处理：{filename}")
