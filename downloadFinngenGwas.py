#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
批量下载 FinnGen R12 GWAS summary statistics（全量 PheWAS 用）

用法：在 Slurm 作业或交互终端里直接
  python downloadFinngenGwas.py
即可把 manifest 里所有 HTTPS 路径下载到 download_directory 目录。

如需并发、断点续传或只下载部分表型，可自行扩展。
"""

import os, requests, pandas as pd

# ======== 1. 路径参数 ========
manifest_path     = "/cluster2/lzhang/finngen_R12_manifest_1000cases.tsv"   # R12 manifest 文件
download_directory = "/cluster2/lzhang/Finngengwas"               # 下载保存目录
os.makedirs(download_directory, exist_ok=True)

# ======== 2. 下载函数 ========
def download_files_from_dataframe(df, url_column, save_directory):
    """逐行下载 df[url_column] 指定的文件到 save_directory"""
    for url in df[url_column]:
        if pd.isna(url):
            continue
        filename = os.path.join(save_directory, os.path.basename(url))
        if os.path.exists(filename):
            print(f"[跳过] 已存在 {filename}")
            continue
        try:
            r = requests.get(url, stream=True, timeout=30)
            r.raise_for_status()
            with open(filename, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"[成功] {filename}")
        except Exception as e:
            print(f"[失败] {url} → {e}")

# ======== 3. 主流程 ========
if __name__ == "__main__":
    df = pd.read_csv(manifest_path, sep="\t")      # 读取 manifest
    url_col = "path_https"                         # manifest 中的 HTTPS 列
    if url_col not in df.columns:
        raise ValueError(f"列 {url_col} 不存在，请检查 manifest")
    download_files_from_dataframe(df, url_col, download_directory)

