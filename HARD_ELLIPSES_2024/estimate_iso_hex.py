#!/usr/bin/env python3
import sys,os
import numpy as np
from operator import itemgetter
filen='_lista_dirs_'
opfn='ordparams.dat'
os.system('ls -d X0_* | sort -t _ -k 2 -n > ' + filen)
with open(filen) as f:
    lines=f.readlines()
if len(sys.argv) > 1:
    psi6thr=float(sys.argv[1])
else:
    psi6thr=0.3
IHEXlst=[]
# fract è la frazione di misure da scartare (intese come equilibratura)
for l in lines:
    os.chdir(l.strip('\n'))
    print('Processing dir '+l.strip('\n')+'\n')
    X0 = l.strip('\n').split('_')[-1]
    with open(opfn) as f:
        opl=f.readlines()
    for ll in opl:
        llst=ll.strip('\n').split()
        phi=float(llst[0])
        psi6=float(llst[1])
        if psi6 > psi6thr:
            phiIHEX=phi
            break
    IHEXlst.append([float(X0),phiIHEX])
    os.chdir('..')
IHEXlst.sort(key=itemgetter(0))
with open('isohex.dat','w') as f:
    for l in IHEXlst:
        f.write(str(l[0])+ ' ' + str(l[1]) + '\n')
os.system('rm ' + filen)
