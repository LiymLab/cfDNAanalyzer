#!/usr/bin/python
"""
generating histograms from TLEN flat files.
"""
import sys
import numpy as np
# Zhou Junpeng 20241105 modified
# np.set_printoptions(suppress=True, precision=200)
#end

if len(sys.argv) < 2:
  print("Usage: %s infile" % sys.argv[0])
  sys.exit(1)
for infile in sys.argv[1:]:
  f = open(infile)
  x=[]
  sum = 0
  for line in f:
    if line[0] != '#':
      fields = line.split()
      x.append(abs(float(fields[0])))
  binedges = range(50, 400, 1)
  vect = np.histogram(x, binedges)      
  vectt = vect[0]
  print( "\t ".join( repr(e) for e in vectt) )
  f.close()
