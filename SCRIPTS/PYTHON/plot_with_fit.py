#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.animation as anim

argv =sys.argv
from scipy.optimize import curve_fit
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
def calcrunavg(vec,lst,lst2):
    c1=1.0
    v0 = vec[0]
    sum=0.0
    for v in vec:
        val=(v-v0)/v0
        sum += val
        lst.append(sum/c1)
        lst2.append(val)
        c1+=1.0
if len(argv) < 2:
    print("Please supply a file with data to plot!")
    quit()
def update(i):
    #f, axarr = plt.subplots()
    plt.clf()
    plt.title(r'orthoterphenyl $\rho\approx 1.050 g/cm^3$ $ T=346 K $ N=2352 $\delta t=5 fs$')
    plt.xlabel('t (ps)')
    for fn in fnlst:
        x, y = np.loadtxt(fn, comments=['#'], usecols=(0,1), unpack=True)
        ra=[]
        yp=[]
        ptd=0
        if df > 0.0:
            ptd = int(len(x)*df)
        if dp > 0:
            ptd=dp
        x = x[ptd:]
        y = y[ptd:]
        if calcra:
            calcrunavg(y, ra, yp)
            y=yp
        sd=calcstddev(y)
        if dofit:
            popt, pcov = curve_fit(func, x, y)
        av=np.mean(y)
        sd=np.std(y)
        plt.ylim(av-sd*sdf,av+sd*sdf)
        plt.plot(x, y, '-', label=fn+' sd=' + str(sd))
        if calcra:
            plt.plot(x, ra,'-r', linewidth=3.0)
        if dofit:
            plt.plot(x, func(x, *popt), '-', label='fit: y=(%G)*x + (%G)' % tuple(popt))
    plt.legend()
    #plt.show()
    plt.draw()

itargs=iter(argv[1:])
sdf=30.0
dp=0
df=0.0
#parse arguments
cc=0
fnlst=[]
keeprunning=False
calcra=False
dofit=False
for a in itargs:
    if a == '-sdf' or a=='--stddevfact':
        sdf=float(next(itargs))
        cc+=1
    elif a == '-df' or a=='--dropfract':
        df=float(next(itargs))
        cc+=1
    elif a == '-dp' or a=='--droppts':
        dp=int(next(itargs))
        cc+=1
    elif a == '-kr' or a=='-keeprun':
        keeprunning=True
    elif a=='-ra' or a == '--runningavg':
        calcra=True
    elif a=='-f' or a=='-fit':
        dofit=True
    elif cc < len(argv)-1:
        fnlst.append(a)
    else:
        print('[ERROR] syntax is the following:')
        print('plot_with_fit.py [-sdf/--stddevfact <this values*stddev will be the yrange> | -df/-dropfract <fraction of points to drop> | -dp/--droppts <number of points to drop>] <data file>')
        quit(1)
    cc=cc+1
    #print('a=', a, ' cc=', cc, ' len(argv)=', len(argv))
if dp > 0 and df > 0:
    print('[ERROR] you have to set either -dp or -df')
    exit(1)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.tick_params(labeltop=True, labelright=True)
if keeprunning:
    a = anim.FuncAnimation(fig, update, repeat=False)
else:
    update(0)
plt.show()
