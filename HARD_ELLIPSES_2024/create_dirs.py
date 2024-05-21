#!/usr/bin/env python3
import os, sys
import numpy as np
fp='lista_press'
if os.path.exists(fp):
    with open(fp) as f:
        l=f.readline()
    T=l.strip('\n').split(' ')
else:
    T=[ str("{:.4f}".format(p)) for p in np.arange(0.5,8.1,0.5) ]
for t in T:
    if t!= '' and (not os.path.exists('P_'+t)):
        os.system('mkdir P_'+t)
