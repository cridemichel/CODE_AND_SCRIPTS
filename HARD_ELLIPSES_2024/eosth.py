#!/usr/bin/env python3
import sys,os
import numpy as np
if len(sys.argv) >= 2:
    X0=float(sys.argv[1])
else:
    print('You must supply an aspect ratio as argument')
    quit()
lst=[]
for nstar in np.arange(0.05, 1.0, 0.01):
    ptonk = nstar/(1.0-nstar)
    psch = 0.5*(ptonk)*(X0-1.0)
    lst.append([nstar,psch])
with open('eosth.dat','w') as f:
    for l in lst:
        f.write(str(l[0])+ ' ' + str(l[1]) + '\n')
