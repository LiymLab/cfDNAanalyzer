import subprocess
import sys

def calculate_coverage(bw_file, region_bed_file, output_file, tool_path, region_name="region"):
    """
    Calculates coverage over specified genomic regions using bigWigAverageOverBed.
    bw_file: Path to the BigWig file.
    region_bed_file: Path to the BED file with specific genomic regions.
    output_file: Path where the output will be saved.
    tool_path: Path to the bigWigAverageOverBed tool.
    region_name: A label for the type of region being analyzed (e.g., TSS, NDR).
    """
    cmd = f"{tool_path} {bw_file} {region_bed_file} {output_file}"
    try:
        subprocess.run(cmd, check=True, shell=True)
        print(f"Coverage calculation for {region_name} completed. Results saved to {output_file}.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {cmd}")
        print(e)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python calculate_coverage.py <bw_file> <region_bed_file> <output_file> <tool_path> <region_name>")
    else:
        _, bw_file, region_bed_file, output_file, tool_path, region_name = sys.argv
        calculate_coverage(bw_file, region_bed_file, output_file, tool_path, region_name)
