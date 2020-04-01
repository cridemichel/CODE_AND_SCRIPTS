#!/usr/bin/env python3
import sys,os
import numpy as np
filen='_lista_dir_'
from operator import itemgetter
os.system('ls -d P_* > ' + filen)
vol=72.3823
N=3360
with open(filen) as f:
    lines=f.readlines()
lst=[]
if len(sys.argv) > 1:
    fract=float(sys.argv[1].strip('\n'))
else:
    fract=0.9
for l in lines:
    os.chdir(l.strip('\n'))
    #press = l.strip('\n').split('_')[-1]
    with open('rho.dat') as f:
        lphi=f.readlines()
    with open('energy.dat') as f:
        lene=f.readlines()
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
    #print('lphi=', lphi)
    nene=len(lene)
    snene=int(fract*nene)
    #print ('press=', press, 'nphi=', nphi, ' snphi=', snphi)
    #print('snphi=', snphi, ' fine=', nphi)
    if snene==nene:
        snene=nene-1
    sublene=lene[snene:nene]
    #print('subl=', subl)
    sumene=0.0
    ccene=0
    for ll in sublene:
        ene=float(ll.strip('\n').split(' ')[1])
        #print('phi=',phi)
        sumene += -ene
        ccene += 1.0
    lst.append([sumphi/float(cc),2.0*sumene/N/float(ccene)])
    #print(sumphi/float(cc),' ',press)
    os.chdir('..')
lst.sort(key=itemgetter(0))
with open('etaf.dat','w') as f:
    for l in lst:
        f.write(str(l[0])+ ' ' + str(l[1]) +'\n')
os.system('rm ' + filen)
