#!/usr/bin/env python

"""Convert BAM files to BigWig file format in a specified region.

Usage:
    bam_to_wiggle.py <BAM file> [<YAML config>]
    [--outfile=<output file name>
     --chrom=<chrom>
     --start=<start>
     --end=<end>
     --normalize]

chrom start and end are optional, in which case they default to everything.
The normalize flag adjusts counts to reads per million.
"""

import os
import sys
import subprocess
import tempfile
from optparse import OptionParser
from contextlib import contextmanager, closing
import yaml
import pysam

def get_program(name, config, ptype="cmd", default=None):
    """Retrieve program information from the configuration."""
    # Directly access the program path from the provided config dictionary
    if 'program' in config and name in config['program']:
        program_path = config['program'][name]
        if ptype == "cmd":
            return program_path  # Assuming the path is directly to the command
        elif ptype == "dir":
            return os.path.dirname(program_path)  # If needing directory, get dirname
    else:
        # Fallback if not specified, use default if provided
        if default is not None:
            return default
        else:
            raise ValueError(f"Program configuration for {name} not found and no default specified.")

def _get_program_cmd(name, config, default):
    """Retrieve commandline of a program."""
    if config is None:
        return name
    elif isinstance(config, str):
        return config
    elif "cmd" in config:
        return config["cmd"]
    elif default is not None:
        return default
    else:
        return name

def _get_program_dir(name, config):
    """Retrieve directory for a program (local installs/java jars)."""
    if config is None:
        raise ValueError("Could not find directory in config for %s" % name)
    elif isinstance(config, str):
        return config
    elif "dir" in config:
        return config["dir"]
    else:
        raise ValueError("Could not find directory in config for %s" % name)

def main(bam_file, config_file=None, chrom='all', start=0, end=None,
         outfile=None, normalize=False, use_tempfile=False):
    if config_file:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
    else:
        config = {"program": {"ucsc_bigwig": "wigToBigWig"}} 

    if outfile is None:
        outfile = f"{os.path.splitext(bam_file)[0]}.bigwig"
    if start > 0:
        start = int(start) - 1
    if end is not None:
        end = int(end)
    regions = [(chrom, start, end)]
    if os.path.abspath(bam_file) == os.path.abspath(outfile):
        print("Bad arguments, input and output files are the same.", file=sys.stderr)
        sys.exit(1)
    if not (os.path.exists(outfile) and os.path.getsize(outfile) > 0):
        if use_tempfile:
            out_handle = tempfile.NamedTemporaryFile(delete=False)
            wig_file = out_handle.name
        else:
            wig_file = f"{os.path.splitext(outfile)[0]}.wig"
            out_handle = open(wig_file, "w")
        with closing(out_handle):
            chr_sizes, wig_valid = write_bam_track(bam_file, regions, config, out_handle, normalize)
        try:
            if wig_valid:
                convert_to_bigwig(wig_file, chr_sizes, config, outfile)
        finally:
            os.remove(wig_file)

@contextmanager
def indexed_bam(bam_file, config):
    if not os.path.exists(f"{bam_file}.bai"):
        pysam.index(bam_file)
    sam_reader = pysam.Samfile(bam_file, "rb")
    yield sam_reader
    sam_reader.close()

def write_bam_track(bam_file, regions, config, out_handle, normalize):
    out_handle.write(f"track {' '.join(['type=wiggle_0', 'name=' + os.path.splitext(os.path.split(bam_file)[-1])[0], 'visibility=full'])}\n")
    normal_scale = 1e6
    is_valid = False
    with indexed_bam(bam_file, config) as work_bam:
        total = sum(1 for r in work_bam.fetch() if not r.is_unmapped) if normalize else None
        sizes = list(zip(work_bam.references, work_bam.lengths))
        if len(regions) == 1 and regions[0][0] == "all":
            regions = [(name, 0, length) for name, length in sizes]
        for chrom, start, end in regions:
            if end is None and chrom in work_bam.references:
                end = work_bam.lengths[work_bam.references.index(chrom)]
            assert end is not None, "Could not find %s in header" % chrom
            out_handle.write(f"variableStep chrom={chrom}\n")
            for col in work_bam.pileup(chrom, start, end):
                n = (float(col.n) / total * normal_scale) if normalize else col.n
                out_handle.write(f"{col.pos + 1} {n:.1f}\n")
                is_valid = True
    return sizes, is_valid

def convert_to_bigwig(wig_file, chr_sizes, config, bw_file=None):
    if not bw_file:
        bw_file = f"{os.path.splitext(wig_file)[0]}.bigwig"
    size_file = f"{os.path.splitext(wig_file)[0]}-sizes.txt"
    with open(size_file, "w") as out_handle:
        for chrom, size in chr_sizes:
            out_handle.write(f"{chrom}\t{size}\n")
    try:
        cl = [get_program("ucsc_bigwig", config, default="wigToBigWig"), wig_file, size_file, bw_file]
        subprocess.check_call(cl)
    finally:
        os.remove(size_file)
    return bw_file

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-o", "--outfile", dest="outfile")
    parser.add_option("-c", "--chrom", dest="chrom")
    parser.add_option("-s", "--start", dest="start")
    parser.add_option("-e", "--end", dest="end")
    parser.add_option("-n", "--normalize", dest="normalize", action="store_true", default=False)
    parser.add_option("-t", "--tempfile", dest="use_tempfile", action="store_true", default=False)
    (options, args) = parser.parse_args()
    if len(args) not in [1, 2]:
        print("Incorrect arguments")
        print(__doc__)
        sys.exit()
    kwargs = dict(
        outfile=options.outfile,
        chrom=options.chrom or 'all',
        start=options.start or 0,
        end=options.end,
        normalize=options.normalize,
        use_tempfile=options.use_tempfile)
    main(*args, **kwargs)
