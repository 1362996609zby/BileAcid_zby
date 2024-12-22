# BileAcid_zby
1. Download verified sequences from UniProt: Search for Bai genes and filter out sequences with Reviewed. Download these sequences as high-quality base sequences.
2. Download DNA sequences of baip, baio, and baij from ((https://doi.org/10.1080/19490976.2022.2132903) and convert them into amino acid sequences using the expasy tool, thereby obtaining all verified amino acid sequence files of Bai.
3. Download unverified sequences of the Bai family from NCBI using the target gene name or annotation keyword search. Compare these sequences with UniProt verified sequences using BLAST and filter out sequences with high similarity (e.g. >60%).
4. Use the tool CD-HIT to remove redundant sequences to avoid duplications causing interference to the model, identity_threshold=0.9.
5. Use the alignment results as input to build the HMM model using hmmbuild.
