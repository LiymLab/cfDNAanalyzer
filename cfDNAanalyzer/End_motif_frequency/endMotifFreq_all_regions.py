import argparse
from Bio import SeqIO
from collections import defaultdict
import math
import itertools
from pybedtools import BedTool
from collections import Counter

# Calculate Motif Diversity Score (MDS)
def calculate_mds(motif_frequencies):
    mds = sum(-p * math.log(p) / math.log(256) for p in motif_frequencies if p > 0)
    return mds

# Find intersection of BED1 and BED2
def intersect_beds(bed1_file, bed2_file):
    bed1 = BedTool(bed1_file)
    bed2 = BedTool(bed2_file)

    bed3 = bed1.intersect(bed2, wa=True, wb=True)
    # print(bed3)
    bed3_regions = []
    for region in str(bed3).strip().split('\n'):
        cols = region.split('\t')
        if len(cols) >= 6:  
            new_region = cols[3:6] + cols[0:3]
            # print(new_region)
            bed3_regions.append(new_region)  # Append the new region to bed3_regions
            
    return bed3_regions
    

# Process BED3 regions to calculate overall motif frequency and MDS
# def process_bed3_and_calculate_mds(bed3_regions, fasta_file):
#     # total_motifs_per_bed1_line = defaultdict(lambda: defaultdict(int))
#     # total_bases_per_bed1_line = defaultdict(int)
#     # fasta_db = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
# 
#     # for region in bed3_regions:
#     #     print("Processing region", region)
#     #     chr_name, start, end, bed1_line = region[0], int(region[1]), int(region[2]), tuple(region[3:])
#     #     seq = str(fasta_db[chr_name].seq[start:end]).upper()
#     # 
#     #     if len(seq) >= 4:
#     #         for i in range(len(seq) - 3):
#     #             subseq = seq[i:i + 4]
#     #             total_motifs_per_bed1_line[bed1_line][subseq] += 1
#     #             total_bases_per_bed1_line[bed1_line] += 1
#     total_motifs_per_bed1_line = defaultdict(Counter)
#     total_bases_per_bed1_line = defaultdict(int)
#     fasta_db = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
#     
#     for region in bed3_regions:
#         print("Processing region", region)
#         chr_name, start, end, bed1_line = region[0], int(region[1]), int(region[2]), tuple(region[3:])
#         seq = str(fasta_db[chr_name].seq[start:end]).upper()
# 
#         if len(seq) >= 4:
#             total_bases_per_bed1_line[bed1_line] += len(seq) - 3  # Increment total bases by the number of possible 4-mers
#             total_motifs_per_bed1_line[bed1_line].update(seq[i:i + 4] for i in range(len(seq) - 3))
# 
#     mds_results = {}
#     frequency_results = []
#     all_motifs = {"".join(motif): 0 for motif in itertools.product("ACGT", repeat=4)}  # All possible 4-mer motifs
def process_bed3_and_calculate_mds(bed3_regions, fasta_file):
    total_motifs_per_bed1_line = defaultdict(Counter)
    total_bases_per_bed1_line = defaultdict(int)
    fasta_db = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    for region in bed3_regions:
        # print("Processing region", region)
        chr_name = region[0]
        start, end = int(region[1]), int(region[2])
        bed1_line = tuple(region[3:])

        seq = str(fasta_db[chr_name].seq[start:end]).upper()
        seq_length = len(seq)

        if seq_length >= 4:
            total_bases_per_bed1_line[bed1_line] += seq_length - 3  # Increment total bases by the number of possible 4-mers
            for i in range(seq_length - 3):
                total_motifs_per_bed1_line[bed1_line][seq[i:i + 4]] += 1

    mds_results = {}
    frequency_results = []
    all_motifs = {''.join(motif): 0 for motif in itertools.product("ACGT", repeat=4)}

    # Calculate MDS and frequency for each BED1 line
    for bed1_line, motifs in total_motifs_per_bed1_line.items():
        # print("Processing BED1 line:", bed1_line)
        total_bases = total_bases_per_bed1_line[bed1_line]

        # Calculate motif frequencies
        motif_frequencies = [count / total_bases for count in motifs.values()] if total_bases > 0 else [0]

        mds = calculate_mds(motif_frequencies)
        mds_results[bed1_line] = mds

        # Fill frequency_results with existing motifs
        for motif, count in motifs.items():
            freq = count / total_bases if total_bases > 0 else 0
            frequency_results.append((bed1_line, motif, freq))

        # Add missing motifs with frequency 0
        for motif in all_motifs.keys():
            if motif not in motifs:
                frequency_results.append((bed1_line, motif, 0))
                
    # print(mds_results)
    # print(frequency_results[:10])
    return mds_results, frequency_results

def main(bed1_file, bed2_file, fasta_file, mds_file, freq_file):
    # Step 1: Find intersection between BED1 and BED2
    bed3_regions = intersect_beds(bed1_file, bed2_file)

    # Step 2: Process BED3 regions to calculate overall motif frequency and MDS
    mds_results, frequency_results = process_bed3_and_calculate_mds(bed3_regions, fasta_file)

    # Step 3: Write MDS results to mds_file in the specified format
    with open(mds_file, 'w') as out_mds_file:
        for bed1_line, mds in mds_results.items():
            chr_name, start, end = bed1_line
            out_mds_file.write(f"{chr_name}\t{start}\t{end}\t{mds:.6f}\n")  # Formatting MDS to 6 decimal places

    # Step 4: Write motif frequency results to freq_file in the specified format
    with open(freq_file, 'w') as out_freq_file:
        for bed1_line, motif, freq in frequency_results:
            chr_name, start, end = bed1_line
            out_freq_file.write(f"{chr_name}\t{start}\t{end}\t{motif}\t{freq:.6f}\n")  # Formatting frequency to 6 decimal places

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Intersect two BED files and calculate overall motif frequency and MDS, output results with BED1 lines.')
    parser.add_argument('-b1', '--bed1', required=True, help='First BED file.')
    parser.add_argument('-b2', '--bed2', required=True, help='Second BED file.')
    parser.add_argument('-f', '--fa', required=True, help='FASTA file of the reference genome.')
    parser.add_argument('-o1', '--output_mds', required=True, help='Output filename for the MDS.')
    parser.add_argument('-o2', '--output_freq', required=True, help='Output filename for the motif frequencies.')
    args = parser.parse_args()
    main(args.bed1, args.bed2, args.fa, args.output_mds, args.output_freq)
