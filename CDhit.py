import pandas as pd
import subprocess
# conda activate cdhit_env #python3 CDhit.py
# 1. 从Excel文件中读取序列
#def read_excel(file_path):
    #"""从Excel文件中读取序列，返回字典格式的序列数据."""
    #df = pd.read_excel(file_path)  # 读取Excel文件
    #sequences = {}

    # 遍历Excel表格，取出gene列和seq列
    #for index, row in df.iterrows():
        #gene = row['gene']
        #seq = row['seq']
        # 只处理非空的gene和seq
        #if pd.notna(gene) and pd.notna(seq):
            #sequences[gene] = seq
    #return sequences

# 2. 创建FASTA文件
#def write_fasta(sequences, output_file):
    #"""写入序列到FASTA文件."""
    #with open(output_file, 'w') as f:
        #for gene_id, sequence in sequences.items():
            #f.write(f">{gene_id}\n")
            #f.write(f"{sequence}\n")
    #print(f"FASTA文件已创建: {output_file}")

# 3. 调用CD-HIT进行去冗余聚类
def run_cd_hit(input_fasta, output_fasta, identity_threshold=0.9):
    """运行CD-HIT，进行去冗余聚类."""
    # CD-HIT命令
    cmd = [
        'cd-hit', 
        '-i', input_fasta,  # 输入FASTA文件
        '-o', output_fasta,  # 输出聚类后的FASTA文件
        '-c', str(identity_threshold)  # 序列相似度阈值
    ]

    # 调用CD-HIT
    try:
        subprocess.run(cmd, check=True)
        print(f"CD-HIT去冗余完成，结果文件: {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"运行CD-HIT时发生错误: {e}")

# 主程序
if __name__ == "__main__":
    # Excel文件路径
    #excel_file = "K07652 ismej.xlsx"
    
    # 从Excel中读取序列
    #sequences = read_excel(excel_file)
    
    # 定义输入输出文件名
    input_fasta = "BaiO_filtered.fasta"
    output_fasta = "result/BaiO_base.fasta"
    
    # 生成FASTA文件
    #write_fasta(sequences, input_fasta)
    
    # 运行CD-HIT进行去冗余
    run_cd_hit(input_fasta, output_fasta, identity_threshold=0.9)
