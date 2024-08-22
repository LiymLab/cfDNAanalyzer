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

def process_sequence(seq, motifs):
    total = 0
    motif_counts = defaultdict(int)
    if len(seq) < 4:
        return motif_counts, total
    
    for i in range(len(seq) - 3):
        motif = seq[i:i+4]
        if motif in motifs:
            motif_counts[motif] += 1
            total += 1
        rev_comp = motif[::-1].translate(str.maketrans("ATGC", "TACG"))
        if rev_comp in motifs:
            motif_counts[rev_comp] += 1
            total += 1
    
    return motif_counts, total

def calculate_region_mds(motif_counts, total):
    motif_frequencies = {motif: count / total if total > 0 else 0 for motif, count in motif_counts.items()}
    return calculate_mds(motif_frequencies.values())

def main(bed_file, fasta_file, output_folder):
    fasta_db = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    motifs = set()
    bases = ["A", "T", "G", "C"]
    for i in bases:
        for j in bases:
            for k in bases:
                for n in bases:
                    motif = i + j + k + n
                    motifs.add(motif)
                    motifs.add(motif[::-1].translate(str.maketrans("ATGC", "TACG")))
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)
    
    # For individual regions
    with open(bed_file, 'r') as bed:
        for idx, line in enumerate(bed):
            if line.startswith('#'): continue
            parts = line.strip().split()
            chr_name, start, end = parts[0], int(parts[1]) + 1, int(parts[2])
            seq = str(fasta_db[chr_name].seq[start-1:end]).upper()
            motif_counts, total = process_sequence(seq, motifs)
            
            # Calculate MDS for this region
            mds = calculate_region_mds(motif_counts, total)
            
            # Output motif frequencies and MDS for this region
            region_output_file = os.path.join(output_folder, f"region_{idx + 1}_motif_frequency_and_mds.txt")
            with open(region_output_file, 'w') as file:
                file.write(f"MDS\t{mds}\n")
                file.write("\n")
                file.write("Motif\tFrequency\n")
                for motif, count in motif_counts.items():
                    file.write(f"{motif}\t{count / total}\n")
    
    # For all regions combined
    motif_counts_all = defaultdict(int)
    total_all = 0
    
    with open(bed_file, 'r') as bed:
        for line in bed:
            if line.startswith('#'): continue
            parts = line.strip().split()
            chr_name, start, end = parts[0], int(parts[1]) + 1, int(parts[2])
            seq = str(fasta_db[chr_name].seq[start-1:end]).upper()
            motif_counts, total = process_sequence(seq, motifs)
            
            # Accumulate motif counts for all regions
            for motif, count in motif_counts.items():
                motif_counts_all[motif] += count
            total_all += total
    
    # Calculate motif frequencies for all motifs
    motif_frequencies_all = {motif: count / total_all if total_all > 0 else 0 for motif, count in motif_counts_all.items()}
    mds_all = calculate_mds(motif_frequencies_all.values())
    
    # Output all motifs frequencies to a single file
    all_motifs_file = os.path.join(output_folder, "all_motifs_frequency.txt")
    with open(all_motifs_file, 'w') as file:
        file.write("Motif\tFrequency\n")
        for motif, freq in motif_frequencies_all.items():
            file.write(f"{motif}\t{freq}\n")
    
    # Output MDS for all motifs
    mds_file = os.path.join(output_folder, "all_motifs_mds.txt")
    with open(mds_file, 'w') as file:
        file.write(f"MDS\t{mds_all}\n")
    
    print(f"Individual region motif frequencies and MDS, and overall motif frequencies with MDS have been outputted to {output_folder}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate the frequency of 4-mer motifs from a BED file and MDS.')
    parser.add_argument('-b', '--bed', required=True, help='BED file containing regions of interest.')
    parser.add_argument('-f', '--fa', required=True, help='FASTA file of the reference genome.')
    parser.add_argument('-o', '--output', default="End_Motif_Freq", help='Output folder for the results.')
    args = parser.parse_args()
    main(args.bed, args.fa, args.output)
