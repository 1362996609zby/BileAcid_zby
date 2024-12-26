import os
import subprocess
from Bio import SeqIO
# conda activate cdhit_env #python3 CDhit_batch.py
def merge_fasta_files(input_dir, merged_fasta):
    """将多个FASTA文件合并成一个大的FASTA文件."""
    with open(merged_fasta, 'w') as output_handle:
        # 遍历输入文件夹中的所有.fasta文件
        for fasta_file in os.listdir(input_dir):
            if fasta_file.endswith(".fasta") or fasta_file.endswith(".fa"):
                fasta_path = os.path.join(input_dir, fasta_file)
                # 将每个FASTA文件的内容追加到目标合并文件
                with open(fasta_path, 'r') as input_handle:
                    output_handle.write(input_handle.read())
    print(f"所有FASTA文件已合并为：{merged_fasta}")

def run_cd_hit(input_fasta, output_fasta, identity_threshold=0.9):
    """运行CD-HIT，进行去冗余聚类."""
    cmd = [
        'cd-hit', 
        '-i', input_fasta,  # 输入FASTA文件
        '-o', output_fasta,  # 输出聚类后的FASTA文件
        '-c', str(identity_threshold)  # 序列相似度阈值
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"CD-HIT去冗余完成，结果文件: {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"运行CD-HIT时发生错误: {e}")

def process_fasta_files(input_dir, output_dir, identity_threshold=0.9):
    """处理多个FASTA文件，合并并运行CD-HIT去冗余."""
    merged_fasta = os.path.join(output_dir, "merged_fasta.fasta")
    
    # 1. 合并多个FASTA文件为一个
    merge_fasta_files(input_dir, merged_fasta)
    
    # 2. 运行CD-HIT进行去冗余
    output_fasta = os.path.join(output_dir, "result_after_cdhit.fasta")
    run_cd_hit(merged_fasta, output_fasta, identity_threshold)

# 主程序
if __name__ == "__main__":
    input_dir = "homologs"  # 存放多个FASTA文件的文件夹
    output_dir = "output"  # 存放结果的文件夹

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 处理FASTA文件
    process_fasta_files(input_dir, output_dir, identity_threshold=0.9)
