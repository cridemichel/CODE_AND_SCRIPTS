#!/usr/bin/env python3
#import matplotlib
#import matplotlib.pyplot as plt
import numpy as np
import sys, os

#matplotlib.rcParams['text.usetex'] = True
argv =sys.argv
absenedr=True
def func(xx, aa, bb):
    return aa * xx  + bb

def calcstddev(vec):
    m=0
    c1=0
    for v in vec:
        m += v
        c1+=1
    m = m/c1
    s=0
    for v in vec:
        s += (v-m)*(v-m)
    s = np.sqrt(s / c1)
    return s

def calcedr(vecx, vecy, lst, xp, lst2,num):
    # calc energy drift according to
    #Yu et al. Chemical Physics 370 (2010) 294–305 (see Fig. 1 therein)
    c1=1.0
    v0 = vecy[0]
    ssum=0.0
    subset=np.arange(0, len(vecy),num)
    print(subset)
    for i in subset:
        vy = vecy[i]
        vx = vecx[i]
        if absenedr==True:
            val=np.abs((vy-v0)/v0)
        else:
            val=(vy-v0)/v0
        ssum += val
        lst.append(ssum/c1)
        lst2.append(val)
        xp.append(vx)
        c1+=1.0

def calcra(vec, lst):
    c1=1.0
    ssum=0.0
    for v in vec:
        ssum += v
        lst.append(ssum/c1)
        c1+=1.0

#if len(argv) < 2:
#    print("Please supply a file with data to plot!")
#    quit()
#parse arguments
ra=[1.0]
runavg=False
dofit=False
calcenedrift=False
if len(sys.argv) == 1:
    print("[ERROR] you have to supply a file with a list of directories")
    quit()
fn=sys.argv[1]
if len(sys.argv) > 2:
    elems=int(sys.argv[2])
else:
    elems=1
with open(fn) as f:
    listadirs=f.readlines()
yp=[1.0]
yp2=[1.0]
xp=[1.0]
for ef in listadirs:
    olddir=os.getcwd()
    os.chdir(ef.strip('\n'))
    x, y = np.loadtxt('energy.dat', comments=['#'], usecols=(0,1), unpack=True)
    yp.clear()
    ra.clear()
    xp.clear()
    calcedr(x, y, ra, xp, yp, elems)
    columns=np.column_stack((xp,yp))
    #print(yp)
    np.savetxt('delE.dat', columns)
    #sd=calcstddev(y)
    #yp2.clear()
    #calcra(y,yp2)
    yp2=ra
    columns=np.column_stack((xp,yp2))
    np.savetxt('runavg_delE.dat', columns)
    os.chdir(olddir)
