#!/bin/bash
#SBATCH --job-name=download_finngen    # 作业名称
#SBATCH --output=download_output.log  # 输出日志文件
#SBATCH -c 2                          # 每个任务使用2个CPU核
#SBATCH --mem=48G                    # 作业总共使用的内存
#SBATCH --partition=superlong         # 使用指定的作业队列
#SBATCH --exclude=cn2                 # 排除 cn2 节点

module load anaconda3/2020-07
module load singularity/3.6.1
source activate Mainwork 

# 执行的命令
python downloadFinngenGwas.py

