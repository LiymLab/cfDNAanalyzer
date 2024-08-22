#!/usr/bin/python
"""
OCF calculations by M.S. Esfahani 2019-08-23 <shahrokh@stanford.edu> Alizadeh and Diehn labs- Stanford University
"""
import sys
import numpy as np
tss = sys.argv[2]
for infile in sys.argv[1:2]:
  f = open(infile)
  up=[]
  down=[]
  sum = 0
  for line in f:
   # up = [0]
   # down = [0]
    if line[0] != '#':
      fields = line.split()
      up.append(float(fields[1])-float(tss))
      down.append(float(fields[1])+abs(float(fields[2]))-1-float(tss))
  up.append(0)
  down.append(0)
  binedges = range(-300,300, 1)
  vect_up = np.histogram(up, binedges)      
  vect_down = np.histogram(down, binedges)
  vectt_up = vect_up[0]
  vectt_down = vect_down[0]
  vv_diff = (vectt_up) - (vectt_down)
  vv_tot = vectt_up + vectt_down
  ocf1 = float(np.sum(-vv_diff[230:250])) + float(np.sum(vv_diff[350:370]))
  ocf2 = float(np.sum(-vv_diff[101:301])) + float(np.sum(vv_diff[301:501]))
  ocfleft = float(np.sum(vv_diff[101:301]))
  ocfleft_total = float(np.sum(vv_tot[101:301]))
  ocfright = float(np.sum(vv_diff[301:501]))
  ocfright_total = float(np.sum(vv_tot[301:501]))
  ocf = [ocf1,ocf2,ocfleft,ocfleft_total,ocfright,ocfright_total]
  if len(ocf) == 0:
    ocf = [0,0,0,0,0,0] # in case the code did not work as expected!!
  print( "\t ".join( repr(e) for e in ocf ))
  #print(*ocf, sep='\t')
  #print( "\t ".join( repr(e) for e in ocf ))
  f.close()
