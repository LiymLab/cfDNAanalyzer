#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *02.03.2015
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict 

import pysam
import gzip
from bx.intervals.intersection import Intersecter, Interval
import random

def isSoftClipped(cigar):
  #Op BAM Description
  #M 0 alignment match (can be a sequence match or mismatch)
  #I 1 insertion to the reference
  #D 2 deletion from the reference
  #N 3 skipped region from the reference
  #S 4 soft clipping (clipped sequences present in SEQ)
  #H 5 hard clipping (clipped sequences NOT present in SEQ)
  #P 6 padding (silent deletion from padded reference)
  #= 7 sequence match
  #X 8 sequence mismatch
  for (op,count) in cigar:
    if op in [4,5,6]: return True
  return False

def aln_length(cigarlist):
    tlength = 0
    for operation, length in cigarlist:
        if operation in [0, 2, 3] or operation >= 6:
            tlength += length
    return tlength

def readIterator(filenames,options):
    if options.pipe or (len(filenames) == 0):
        if options.verbose: sys.stderr.write("Reading from pipe...\n")
        input_file = pysam.AlignmentFile("-", "r")
        for read in input_file:
            yield read
        input_file.close()
    else:
        for bamfile in filenames:
            if options.verbose: sys.stderr.write("Reading %s\n" % bamfile)
            if os.path.exists(bamfile) and (os.path.exists(bamfile + ".bai") or os.path.exists(bamfile.replace(".bam", ".bai"))):
                input_file = pysam.AlignmentFile(bamfile, "rb")
                if chrom and start and end:
                    for read in input_file.fetch(chrom, start - 1, end):
                        yield read
                else:
                    for read in input_file:
                        yield read
                input_file.close()
                    

def write_output(file_path, data, start, end, chromosome, data_extractor):
    header = f"fixedStep chrom=chr{chromosome} start={start} step=1\n"
    mode = 'wt' if file_path in ["", "OFF", "stdout"] else 'wb'
    output = sys.stdout if file_path in ["", "OFF", "stdout"] else gzip.open(file_path, mode)
    
    if mode == 'wb':
        output.write(header.encode())
        for pos in range(start, end + 1):
            output.write(f"{data_extractor(data, pos)}\n".encode())
    else:
        output.write(header)
        for pos in range(start, end + 1):
            output.write(f"{data_extractor(data, pos)}\n")

    if file_path not in ["", "OFF", "stdout"]:
        output.close()

parser = OptionParser()

parser.add_option("-r","--region", dest="region", help="Region to be looked up (def chr12:34,443,233-34,453,733)",default="chr12:34,443,233-34,453,733")
parser.add_option("-p","--pipe", dest="pipe", help="Do not fetch, stream data (def Off)",default=False,action="store_true")
parser.add_option("-l","--lengthSR", dest="lengthSR", help="Length of full reads (default 76)",default=76,type="int")
parser.add_option("-m","--merged", dest="merged", help="Assume reads are merged (default Off)",default=False,action="store_true")
parser.add_option("-t","--trimmed", dest="trimmed", help="Assume reads are trimmed (default Off)",default=False,action="store_true")
parser.add_option("-w","--protection", dest="protection", help="Base pair protection assumed for elements (default 120)",default=120,type="int")
parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
parser.add_option("-c","--coverage", dest="coverage", help="Coverage output stream (default 'OFF', 'OFF' for no output)",default='OFF')
parser.add_option("-n","--metric", dest="metric", help="Metric output stream (default stdout, 'OFF' for no output)",default="")
parser.add_option("-s","--starts", dest="starts", help="Read starts output stream (default 'OFF', 'OFF' for no output)",default='OFF')
parser.add_option("--random", dest="random", help="Add +/- 5bp wiggle to read ends to obscure fragmentation patterns (default OFF)",default=False,action="store_true")
parser.add_option("--downsample", dest="downsample", help="Ratio to down sample reads (default OFF)",default=None,type="float")
(options, args) = parser.parse_args()

chrom,start,end = None,None,None
outchrom = None
try:
    chrom = options.region.split(':')[0]
    start, end = list(map(int, options.region.split(':')[1].replace(",", "").split("-")))
    outchrom = chrom.replace("chr", "")
    if options.verbose:
        sys.stderr.write("Coordinates parsed: Chrom %s Start %d End %d\n" % (chrom, start, end))
except Exception as e:
    sys.stderr.write("Error parsing region: %s\n" % str(e))
    sys.exit(1)

posRange = defaultdict(lambda:[0,0,0])
filteredReads = Intersecter()
for read in readIterator(args,options):
  if read.is_duplicate or read.is_qcfail or read.is_unmapped: continue
  if isSoftClipped(read.cigar): continue

    
  
  if read.is_paired:
    if read.mate_is_unmapped: continue
    if read.rnext != read.tid: continue
    if read.is_read1 or (not options.pipe and read.is_read2 and read.pnext+read.qlen < start):
      if read.isize == 0: continue
      if options.downsample != None and random.random() >= options.downsample: continue
      if options.random:
        rstart = min(read.pos,read.pnext)+1+random.randint(-5,5) # 1-based
        rend = rstart+abs(read.isize)-1+random.randint(-5,5) # end included
      else:
        rstart = min(read.pos,read.pnext)+1 # 1-based
        rend = rstart+abs(read.isize)-1 # end included
      
      filteredReads.add_interval(Interval(rstart,rend))
      #print read.qname,rstart,rend,rend-rstart,abs(read.isize)
      for i in range(rstart,rend+1):
        if i >= start and i <= end:
          posRange[i][0]+=1
      if rstart >= start and rstart <= end:
        posRange[rstart][1]+=1
      if rend >= start and rend <= end:
        posRange[rend][1]+=1
  else:
    if options.downsample != None and random.random() >= options.downsample: continue
    if options.random:
      rstart = read.pos+1+random.randint(-5,5) # 1-based
      rend = rstart+aln_length(read.cigar)-1+random.randint(-5,5) # end included
    else:
      rstart = read.pos+1 # 1-based
      rend = rstart+aln_length(read.cigar)-1 # end included
    
    filteredReads.add_interval(Interval(rstart,rend))
    #print read.qname,rstart,rend,rend-rstart,aln_length(read.cigar)
    for i in range(rstart,rend+1):
      if i >= start and i <= end:
        posRange[i][0]+=1
    if ((options.merged or read.qname.startswith('M_')) or ((options.trimmed or read.qname.startswith('T_')) and read.qlen <= options.lengthSR-10)):
      if (rstart >= start and rstart <= end):
        posRange[rstart][1]+=1
      if rend >= start and rend <= end:
        posRange[rend][1]+=1
    elif read.is_reverse:
      if rend >= start and rend <= end:
        posRange[rend][1]+=1
    else:
      if (rstart >= start and rstart <= end):
        posRange[rstart][1]+=1

protection = options.protection//2
for pos in range(start,end+1):
  rstart,rend = pos-protection,pos+protection
  gcount,bcount = 0,0
  for read in filteredReads.find(rstart,rend):
     if (read.start > rstart) or (read.end < rend): bcount +=1
     else: gcount +=1
  posRange[pos][2]+=gcount-bcount

write_output(options.coverage, posRange, start, end, outchrom, lambda data, pos: data[pos][0])
write_output(options.starts, posRange, start, end, outchrom, lambda data, pos: data[pos][1])
write_output(options.metric, posRange, start, end, outchrom, lambda data, pos: data[pos][2])
