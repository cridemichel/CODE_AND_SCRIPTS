#!/usr/bin/env python3
import os, sys
import numpy as np
fp='lista_press'
if len(sys.argv) >= 2:
    minp=sys.argv[1]
else:
    minp=4.9
if len(sys.argv) >= 3:
    maxp=sys.argv[2]
else:
    maxp=15.1
if len(sys.argv) >= 4:
    dp=sys.argv[3]
else:
    dp=0.5
if os.path.exists(fp):
    with open(fp) as f:
        l=f.readline()
    T=l.strip('\n').split(' ')
else:
    T=[ str("{:.4f}".format(p)) for p in np.arange(float(minp),float(maxp),float(dp)) ]
for t in T:
    if t!= '' and (not os.path.exists('P_'+t)):
        os.system('mkdir P_'+t)
