#!/usr/bin/env python3
import sys,os
import numpy as np
filen='_lista_dir_'
from operator import itemgetter
os.system('ls -d P_*[^a-z] | sort -t _ -k 2 -n > ' + filen)
vol=1.0
pfact=2.0
nfact=2.0
with open(filen) as f:
    lines=f.readlines()
lst=[]
if len(sys.argv) > 1:
    fract=float(sys.argv[1].strip('\n'))
else:
    fract=0.8
for l in lines:
    os.chdir(l.strip('\n'))
    press = l.strip('\n').split('_')[-1]
    with open('rho.dat') as f:
        lphi=f.readlines()
    #print('lphi=', lphi)
    nphi=len(lphi)
    snphi=int(fract*nphi)
    #print ('press=', press, 'nphi=', nphi, ' snphi=', snphi)
    #print('snphi=', snphi, ' fine=', nphi)
    if snphi==nphi:
        snphi=nphi-1
    subl=lphi[snphi:nphi]
    #print('subl=', subl)
    sumphi=0.0
    cc=0
    for ll in subl:
        phi=vol*float(ll.strip('\n').split(' ')[1])
        #print('phi=',phi)
        sumphi+=phi
        cc += 1.0
    nstar = nfact*sumphi/float(cc)
    # pressure of X=1 1D fluid (Tonk fluid)
    ptonk = nstar/(1.0-nstar)
    lst.append([nstar,2.0*float(press)/ptonk-1.0])
    #print(sumphi/float(cc),' ',press)
    os.chdir('..')
lst.sort(key=itemgetter(0))
with open('eos.dat','w') as f:
    for l in lst:
        f.write(str(l[0])+ ' ' + str(l[1]) +'\n')
os.system('rm ' + filen)
