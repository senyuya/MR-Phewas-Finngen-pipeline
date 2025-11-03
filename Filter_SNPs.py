import os
import pandas as pd
import gzip

# 默认已有一个snplist.txt文件，为一列snps的纯文本文件（纯rs列表）

# 路径设置
gwas_dir = "/cluster2/lzhang/Finngengwas"
snp_file = "/cluster2/lzhang/snplist.txt"
output_dir = "/cluster2/lzhang/extracted_snps"

# 创建输出文件夹
os.makedirs(output_dir, exist_ok=True)

# 读取SNP列表
with open(snp_file, 'r') as f:
    target_snps = set(line.strip() for line in f if line.strip())

print(f"需要提取的SNP数量: {len(target_snps)}")

# 获取所有GWAS文件
gwas_files = [f for f in os.listdir(gwas_dir) if f.endswith('.gz')]
print(f"找到GWAS文件数量: {len(gwas_files)}")

# 处理每个文件
for i, gwas_file in enumerate(gwas_files):
    try:
        print(f"处理 {i+1}/{len(gwas_files)}: {gwas_file}")
        
        # 读取文件
        filepath = os.path.join(gwas_dir, gwas_file)
        with gzip.open(filepath, 'rt') as f:
            gwas_data = pd.read_csv(f, sep='\t')
        
        # 提取匹配的SNP
        matches = gwas_data[gwas_data['rsids'].isin(target_snps)]
        
        # 如果有匹配的SNP，保存结果
        if not matches.empty:
            output_file = os.path.join(output_dir, f"extracted_{gwas_file.replace('.gz', '.txt')}")
            matches.to_csv(output_file, sep='\t', index=False)
            print(f"  找到 {len(matches)} 个匹配SNP，已保存")
        else:
            print(f"  没有找到匹配的SNP")
            
    except Exception as e:
        print(f"  处理文件出错: {e}")

print("处理完成！")