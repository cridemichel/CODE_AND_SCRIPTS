#!/usr/bin/env python3
import sys,os
import numpy as np
def avgval(fn, fract):
    with open(fn) as f:
        lst=f.readlines()
    #print('lphi=', lphi)
    nlst=len(lst)
    snlst=int(fract*nlst)
    if snlst==nlst:
        snlst=nlst-1
    #print('subl=', subl)
    subl=lst[snlst:nlst]
    sum=0.0
    cc=0
    for ll in subl:
        val=float(ll.strip('\n').split(' ')[1])
        #print('phi=',phi)
        sum+=val
        cc += 1.0
    avgval = sum/float(cc)
    print('fn=', fn, ' avgval= ', avgval, ' cc=', cc)
    return avgval
filen='_lista_dir_'
from operator import itemgetter
if !os.path.exists('lista_dirs'):
    os.system('ls -d P_*[^a-z] | sort -t _ -k 2 -n > ' + filen)
else:
    os.system('cp lista_dirs ' + filen)
pfact=1.0
nfact=1.0
with open(filen) as f:
    lines=f.readlines()
lst=[]
# fract è la frazione di misure da scartare (intese come equilibratura)
if len(sys.argv) > 1:
    fract=float(sys.argv[1].strip('\n'))
else:
    fract=0.8
if len(sys.argv) > 2:
    X0=float(sys.argv[2])
else:
    X0=1.0
vol=X0*3.141592653589793
for l in lines:
    os.chdir(l.strip('\n'))
    print('Processing dir '+l.strip('\n')+'\n')
    press = l.strip('\n').split('_')[-1]
    avgphi=avgval('phi.dat',fract)
    # pressure of X=1 1D fluid (Tonk fluid)
    #ptonk = nstar/(1.0-nstar)
    # Schilling prediction
    #psch = 0.5*(ptonk)*(X0-1.0)
    os.system('ls cnf-* | sort -t - -k 2 -n > _listacnf_')
    os.system('wc -l _listacnf_ | gawk \'{print($1)}\' > _tmp_')
    NF=open('_tmp_', 'r').read().strip('\n')
    print('NF='+NF)
    os.system('rm _tmp_')
    os.system('cat _listacnf_ | tail -n ' + str(int((1.0-fract)*int(NF))) + ' > _sublistacnf_')
    os.system('../../calc_order_params _sublistacnf_ > /dev/null')
    #os.system('rm _listacnf_ _sublistacnf_')
    avgpsi6=avgval('psi6.dat',1.0)
    avgpsiT=avgval('psiT.dat',1.0)
    avgS=avgval('S.dat',1.0)
    lst.append([avgphi,float(avgpsi6),float(avgpsiT),float(avgS)])
    #print(sumphi/float(cc),' ',press)
    os.chdir('..')
lst.sort(key=itemgetter(0))
with open('ordparams.dat','w') as f:
    for l in lst:
        f.write(str(l[0])+ ' ' + str(l[1]) + ' ' + str(l[2]) + ' ' + str(l[3]) +  '\n')
os.system('rm ' + filen)
