#!/usr/bin/env python3
import sys,os
import numpy as np
filen='_lista_dir_'
#from operator import itemgetter
def optd(partsL):
    d=1.05
    dend=1.2
    smop_o=smop_oo=-1.0
    deld=0.005
    fine=0
    while not fine:
        #print('lphi=', lphi)
        #print('snphi=', snphi, ' fine=', nphi)
        sumca=0.0
        sumsa=0.0
        cc=0
        for p in partsL[1:]:
            ap=p.strip('\n').split()
            y=float(ap[1])
            a = 2.0*np.pi*y/d
            sumca+=np.cos(a)
            sumsa+=np.sin(a)
            cc += 1.0
        smop = sumca*sumca/cc/cc + sumsa*sumsa/cc/cc
        #       print('smop_oo=', smop_oo, ' smop_o=', smop_o, ' smop=', smop)
        if smop_oo > 0 and smop_o > 0:
            if smop_o > smop_oo and smop_o > smop:
                #it is maximum
                return d
        smop_oo=smop_o
        smop_o=smop
        #print(d,' ', smop) 
        d=d+deld
        if d >= dend:
            return -1
            fine=1

if len(sys.argv)==3:
    os.system('cp '+sys.argv[2]+' '+filen)
else:
    os.system('ls -d P_*[^a-z] | sort -t _ -k 2 -n > ' + filen)
vol=1.0
with open(filen) as f:
    lines=f.readlines()
#print(lines)
lst=[]
dsm=1.13
if len(sys.argv) > 1:
    fract=float(sys.argv[1].strip('\n'))
else:
    fract=0.9
for l in lines:
    print('processing dir=', l)
    os.chdir(l.strip('\n'))
    press = l.strip('\n').split('_')[-1]
    if fract < 0:
        stri='ls cnf-*[0-9] | sort -t - -k 2 -n | tail -n ' + str(abs(int(fract))) + ' > listacnf'
        #print(stri)
        os.system(stri)
    else:
        os.system('ls cnf-*[0-9] | sort -t - -k 2 -n > listacnf') #crea la lista dei file cnf da utilizzare
    with open('listacnf') as f:  # <==== apre il file con lista dei file
        lf=f.readlines()
    #print('lphi=', lphi)
    nc=len(lf)
    if fract < 0:
        snc=0
    else:
        snc=int(fract*nc)
    #print ('press=', press, 'nphi=', nphi, ' snphi=', snphi)
    #print('snphi=', snphi, ' fine=', nphi)
    if snc==nc:
        snc=nc-1
    subl=lf[snc:nc]
    #print('subl=', subl)
    smop=0.0
    cc2=0.0
    for ll in subl:
        with open(ll.strip('\n')) as acnf: # <===== apre un file cnf
            parts=acnf.readlines()
        dsm=optd(parts)
        if dsm==-1:
            dsm=1.1
        print('dsm=', dsm)
        cc=0.0
        sumca=0.0
        sumsa=0.0
        for p in parts[1:]:
            ap=p.strip('\n').split()
            y=float(ap[1])
            a = 2.0*np.pi*y/dsm
            sumca+=np.cos(a)
            sumsa+=np.sin(a)
            cc += 1.0
        smop += sumca*sumca/cc/cc + sumsa*sumsa/cc/cc
        cc2 += 1.0
#print('phi=',phi)
    lst.append([press,str(smop/cc2)])
    #print(sumphi/float(cc),' ',press)
    os.chdir('..')
#lst.sort(key=itemgetter(0))
with open('smordpar.dat','w') as f:
    for ll in lst:
        f.write(str(ll[0])+ ' ' + ll[1]+'\n')
#os.system('rm ' + filen)
