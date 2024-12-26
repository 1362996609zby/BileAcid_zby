from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
import os
#conda activate cdhit_env #python blast_local_filter.py
def run_blast_local(query_fasta, subject_fasta, output_dir, identity_threshold=60):
    """
    Run BLAST locally to compare unverified sequences against verified sequences.

    Parameters:
    - query_fasta (str): Path to the query FASTA file (e.g., NCBI unverified sequences).
    - subject_fasta (str): Path to the subject FASTA file (e.g., UniProt verified sequences).
    - output_dir (str): Directory to store all output files.
    - identity_threshold (float): Minimum identity percentage to retain.

    Output:
    - Writes filtered sequences and BLAST results to specified directory.
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Step 1: Create BLAST database from the subject FASTA (UniProt verified sequences)
    print("Creating BLAST database from subject FASTA...")
    os.system(f"makeblastdb -in {subject_fasta} -dbtype prot -out {os.path.join(output_dir, 'subject_db')}")

    # Step 2: Run BLASTP with the query FASTA against the subject database
    print("Running BLAST...")
    blast_output_file = os.path.join(output_dir, "blast_results.txt")
    blastp_cline = NcbiblastpCommandline(
        query=query_fasta,
        db=os.path.join(output_dir, "subject_db"),
        out=blast_output_file,
        evalue=0.001,
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    )
    stdout, stderr = blastp_cline()
    print("BLAST completed. Results saved to:", blast_output_file)

    # Step 3: Filter results based on identity threshold
    filtered_sequences = []
    with open(blast_output_file, "r") as blast_results:
        for line in blast_results:
            parts = line.strip().split("\t")
            query_id, subject_id, identity = parts[0], parts[1], float(parts[2])

            if identity >= identity_threshold:
                print(f"Match: {query_id} -> {subject_id} with {identity}% identity")
                filtered_sequences.append(query_id)

    # Step 4: Write filtered sequences to a new FASTA file
    filtered_fasta = os.path.join(output_dir, "filtered_sequences.fasta")
    with open(filtered_fasta, "w") as out_fasta:
        for record in SeqIO.parse(query_fasta, "fasta"):
            if record.id in filtered_sequences:
                SeqIO.write(record, out_fasta, "fasta")

        # Append UniProt verified sequences to the filtered FASTA file
        for record in SeqIO.parse(subject_fasta, "fasta"):
            SeqIO.write(record, out_fasta, "fasta")

    print(f"Filtered sequences (including UniProt verified sequences) saved to {filtered_fasta}")


# Example usage
run_blast_local("ncbi_all/BaiO_ncbi.fasta", "BaiO.fasta", "BaiO_filtered", identity_threshold=60)
