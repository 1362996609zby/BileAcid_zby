import os
import pandas as pd
from Bio import SeqIO
# conda activate hmm_env
# export PATH=/root/miniconda3/envs/hmm_env/bin:$PATH
# python extract_homologs_from_MAG.py
def extract_homologs_from_MAG(excel_file, mag_dir, output_dir):
    """
    Extract homologous sequences from MAG files based on an Excel table with target names and corresponding MAG file names.
    
    Parameters:
    - excel_file (str): Path to the Excel file containing target names, query names (bai gene names), and corresponding MAG file names.
    - mag_dir (str): Directory containing translated MAG sequence files (with .fasta extension).
    - output_dir (str): Directory to save the extracted homologous sequences in .fasta format.
    """
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Step 1: Read the Excel file and extract target names, query names (bai gene names), and MAG file names
    df = pd.read_excel(excel_file)
    target_to_file_and_query = {}

    # Populate a dictionary with target_name, query_name (bai gene) and corresponding MAG file names
    for _, row in df.iterrows():
        target_name = row['target_name']
        bai_gene_name = row['query_name'].replace("_base_alignment", "")  # Remove "_base_alignment" from query_name
        mag_file = row['file'].replace(".txt", "_translated.fasta")  # Modify file extension
        if target_name not in target_to_file_and_query:
            target_to_file_and_query[target_name] = []
        target_to_file_and_query[target_name].append((bai_gene_name, mag_file))

    # Step 2: For each target, extract the homologous sequence from the corresponding MAG file
    for target_name, bai_genes in target_to_file_and_query.items():
        for bai_gene_name, mag_file in bai_genes:
            # Construct the path to the translated MAG file
            mag_path = os.path.join(mag_dir, mag_file)

            if os.path.exists(mag_path):  # Ensure the MAG file exists
                # Step 3: Construct the output file name by adding bai gene name and MAG file name
                output_file = os.path.join(output_dir, f"{target_name}_{bai_gene_name}_{mag_file}_homolog.fasta")

                # Step 4: Read the MAG file and extract the homologous sequence
                with open(mag_path, "r") as mag_fasta, open(output_file, "w") as out_fasta:
                    records = SeqIO.parse(mag_fasta, "fasta")
                    for record in records:
                        if record.id == target_name:  # Match the target sequence id
                            SeqIO.write(record, out_fasta, "fasta")
                print(f"Extracted homologous sequence for {target_name} ({bai_gene_name}) saved to {output_file}")
            else:
                print(f"Warning: MAG file {mag_file} not found in {mag_dir}. Skipping.")

# Example usage
excel_file = "ismej_dis.xlsx"  # Path to your Excel file containing target_name, query_name, and file columns
mag_dir = "MAG_translated_sequences"   # Directory containing the translated MAG files
output_dir = "homologs"  # Directory to save the extracted homologous sequences

extract_homologs_from_MAG(excel_file, mag_dir, output_dir)
