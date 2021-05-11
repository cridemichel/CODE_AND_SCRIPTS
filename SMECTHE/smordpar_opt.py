#!/usr/bin/env python3
import sys,os
import numpy as np
from operator import itemgetter
vol=1.0
cn = sys.argv[1]
#print(lines)
lst=[]
d=0.9
dend=1.4
fine=0
while not fine:
    #print('lphi=', lphi)
    #print('snphi=', snphi, ' fine=', nphi)
    smop=0.0 
    sumca=0.0
    sumsa=0.0
    cc=0
    with open(cn) as acnf:
        parts=acnf.readlines()
    deld=0.005
    dmax=1.5
    for p in parts[1:]:
        ap=p.strip('\n').split()
        y=float(ap[1])
        a = 2.0*np.pi*y/d
        sumca+=np.cos(a)
        sumsa+=np.sin(a)
        cc += 1.0
    smop = sumca*sumca/cc/cc + sumsa*sumsa/cc/cc
    print(d,' ', smop) 
    d=d+deld
    if d >= dend: 
        fine=1
