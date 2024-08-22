# To use: python endMotifFreq.py -b example.bed -f ref.fasta -o output_dir
# conda install biopython
import argparse
from Bio import SeqIO
from collections import defaultdict
import math
import os

# Calculate Motif Diversity Score
def calculate_mds(motif_frequencies):
    mds = sum(-p * math.log(p) / math.log(256) for p in motif_frequencies if p > 0)
    return mds

def main(bed_file, fasta_file, output_folder):
    motifs = defaultdict(int)
    total = 0
    fasta_db = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    bases = ["A", "T", "G", "C"]
    for i in bases:
        for j in bases:
            for k in bases:
                for n in bases:
                    motif = i + j + k + n
                    motifs[motif] = 0
    
    with open(bed_file, 'r') as bed:
        for line in bed:
            if line.startswith('#'): continue
            parts = line.strip().split()
            chr_name, start, end = parts[0], int(parts[1]) + 1, int(parts[2])
            seq = str(fasta_db[chr_name].seq[start-1:end]).upper()
            if seq[:4] in motifs:
                motifs[seq[:4]] += 1
                total += 1
            rev_comp = seq[-4:][::-1].translate(str.maketrans("ATGC", "TACG"))
            if rev_comp in motifs:
                motifs[rev_comp] += 1
                total += 1

    # Calculate motif frequencies and MDS
    motif_frequencies = []
    for motif in motifs:
        freq = motifs[motif] / total if total > 0 else 0
        motif_frequencies.append(freq)
    mds = calculate_mds(motif_frequencies)
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)

    freq_file = os.path.join(output_folder, "motifs_frequency.txt")
    with open(freq_file, 'w') as file:
        for motif, freq in zip(motifs.keys(), motif_frequencies):
            file.write(f"{motif}\t{freq}\n")
    
    mds_file = os.path.join(output_folder, "motifs_mds.txt")
    with open(mds_file, 'w') as file:
        file.write(f"MDS\t{mds}\n")

    print(f"Motif frequencies and MDS have been outputted to {output_folder}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate the frequency of 4-mer motifs from a BED file and MDS.')
    parser.add_argument('-b', '--bed', required=True, help='BED file containing regions of interest.')
    parser.add_argument('-f', '--fa', required=True, help='FASTA file of the reference genome.')
    parser.add_argument('-o', '--output', default="End_Motif_Freq", help='Output folder for the results.')
    args = parser.parse_args()
    main(args.bed, args.fa, args.output)
