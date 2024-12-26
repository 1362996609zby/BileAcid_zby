import os
from Bio.Align.Applications import MafftCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# conda create -n hmm_env python=3.9 -y
# conda activate hmm_env
# export PATH=/root/miniconda3/envs/hmm_env/bin:$PATH
# python Bai_Gene_Hmm.py
def translate_fasta(input_dir, output_dir):
    """
    Translate nucleotide sequences in .fa files to protein sequences.
    
    Parameters:
    - input_dir (str): Directory containing nucleotide .fa files.
    - output_dir (str): Directory to save translated protein sequences as .fasta.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file in os.listdir(input_dir):
        if file.endswith(".fna"):
            input_path = os.path.join(input_dir, file)
            output_path = os.path.join(output_dir, f"{os.path.splitext(file)[0]}_translated.fasta")
            with open(output_path, "w") as output_file:
                for record in SeqIO.parse(input_path, "fasta"):
                    # Translate nucleotide sequence into protein sequence (six frames)
                    for frame in range(3):
                        translated_seq = record.seq[frame:].translate(to_stop=True)
                        translated_record = SeqRecord(
                            translated_seq,
                            id=f"{record.id}_frame{frame + 1}",
                            description=f"Translated frame {frame + 1}"
                        )
                        SeqIO.write(translated_record, output_file, "fasta")
            print(f"Translated {file} to {output_path}")

def construct_and_search_per_bai(bai_sequences_dir, mag_dir, output_dir, mag_translated_dir):
    """
    Construct individual HMM models for each bai gene and search in MAGs.

    Parameters:
    - bai_sequences_dir (str): Directory containing bai gene FASTA files.
    - mag_dir (str): Directory containing MAG FASTA files.
    - output_dir (str): Directory to save the HMM search results.
    - mag_translated_dir (str): Directory containing translated MAG sequences.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for fasta_file in os.listdir(bai_sequences_dir):
        if fasta_file.endswith(".fasta"):
            bai_name = os.path.splitext(fasta_file)[0]
            alignment_output = os.path.join(output_dir, f"{bai_name}_alignment.sto")
            hmm_output = os.path.join(output_dir, f"{bai_name}.hmm")
    # Step 0: Translate MAG sequences if they are in nucleotide format
    translate_fasta(mag_dir, mag_translated_dir)

    for fasta_file in os.listdir(bai_sequences_dir):
        if fasta_file.endswith(".fasta"):
            bai_name = os.path.splitext(fasta_file)[0]
            alignment_output = os.path.join(output_dir, f"{bai_name}_alignment.sto")
            hmm_output = os.path.join(output_dir, f"{bai_name}.hmm")
            # Step 1: Perform multiple sequence alignment using MAFFT
            print(f"Aligning sequences for {bai_name}...")
            mafft_cline = MafftCommandline(input=os.path.join(bai_sequences_dir, fasta_file))
            stdout, stderr = mafft_cline()
            with open(alignment_output, "w") as align_file:
                align_file.write(stdout)
            print(f"Alignment for {bai_name} saved to {alignment_output}.")

            # Step 2: Build HMM model using hmmbuild
            print(f"Building HMM model for {bai_name}...")
            os.system(f"hmmbuild {hmm_output} {alignment_output}")
            print(f"HMM model for {bai_name} saved to {hmm_output}.")

            # Step 3: Search HMM model in MAGs
            bai_output_dir = os.path.join(output_dir, bai_name)
            if not os.path.exists(bai_output_dir):
                os.makedirs(bai_output_dir)

            for mag_file in os.listdir(mag_translated_dir):
                if mag_file.endswith(".fasta"):
                    mag_path = os.path.join(mag_translated_dir, mag_file)
                    output_file = os.path.join(bai_output_dir, f"{os.path.splitext(mag_file)[0]}_hmmsearch.txt")
                    print(f"Searching {bai_name} in {mag_file}...")
                    os.system(f"hmmsearch --tblout {output_file} -E 1e-10 {hmm_output} {mag_path}")
                    print(f"HMM search results for {bai_name} in {mag_file} saved to {output_file}.")

# Example usage
construct_and_search_per_bai("bai_gene_sequences", "iso_MAG_sequences", "HMM_search_results_per_bai", "MAG_translated_sequences")
