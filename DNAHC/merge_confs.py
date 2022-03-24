#!/usr/local/bin/python3
#import numpy as np
import sys
#Usage calc_folding.py <threshold_in_degree> 
fn1=sys.argv[1]
fn2=sys.argv[2]
with open(fn1,encoding='utf-8') as f:
    lines1=f.readlines()
with open(fn2,encoding='utf-8') as f:
    lines2=f.readlines()
#  if len(sys.argv) > 2:
#      thr=float(sys.argv[2])
#  else:
#      thr=15.0
dely1=1.4
dely2=2.1
hdr1=lines1[0].strip('\n').split()
hdr2=lines2[0].strip('\n').split()
N1=int(hdr1[0])
N2=int(hdr2[0])
box1=[float(i) for i in hdr1[1:]]
box2=[float(i) for i in hdr2[1:]]
print(box1,box2)
lines1f=[]
for ll in lines1[1:]:
    lines1f.append([float(i) for i in ll.strip('\n').split()])
lines2f=[]
for ll in lines2[1:]:
    lines2f.append([float(i) for i in ll.strip('\n').split()])
Lx=max(box1[0],box2[0])
Ly=box1[1]+box2[1]+2.0*dely2
Lz=max(box1[2],box2[2])
newconf=[]
Ntot=0
maxy=0.0
miny=0.0
for l in lines1f:
    l[1] = l[1]-box1[1]/2.0-dely1
    if l[1] > -Ly/4.0:
        Ntot+=1
        if l[1] < miny:
            miny=l[1]
        newconf.append(l)
for l in lines2f:
    l[1] = l[1]+box2[1]/2.0+dely1
    if l[1] < Ly/4.0:
        Ntot += 1
        if l[1] > maxy:
            maxy = l[1]
        newconf.append(l)
print('maxy=', maxy, ' miny=', miny, ' maxy-miny=', maxy-miny)
newhdrf=[Ntot, Lx, (maxy-miny)+dely2, Lz]
newhdr=' '.join([str(i) for i in newhdrf ])
#print(newconf)
#quit()
with open("cnf-merged","w",encoding='utf-8') as f:
    f.write(newhdr+'\n')
    for lll in newconf:
        f.write(' '.join([str(i) for i in lll ])+'\n')
