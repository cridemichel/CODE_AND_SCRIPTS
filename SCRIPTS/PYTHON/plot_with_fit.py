#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
argv =sys.argv
from scipy.optimize import curve_fit
def func(x, a, b,):
    return a * x  + b
if len(argv) < 2:
    print("Please supply a file with data to plot!")
    quit()
itargs=iter(argv[1:])
sdf=30.0
dp=0
df=0.0
#parse arguments
cc=0

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
    elif cc==len(argv)-2:
        fn=a
    else:
        print('[ERROR] syntax is the following:')
        print('plot_with_fit.py [-sdf/--stddevfact <this values*stddev will be the yrange> | -df/-dropfract <fraction of points to drop> | -dp/--droppts <number of points to drop>] <data file>')
        quit(1)
    cc=cc+1
    #print('a=', a, ' cc=', cc, ' len(argv)=', len(argv))
if dp > 0 and df > 0:
    print('[ERROR] you have to set either -dp or -df')
    exit(1)
f, axarr = plt.subplots()
axarr.tick_params(labeltop=True, labelright=True)
plt.title(r'orthoterphenyl $\rho\approx 1.106 g/cm^3$ $ T=281 K $ N=2352 $\delta t=5 fs$')
plt.xlabel('t (ps)')
plt.ylabel('energy')
x, y = np.loadtxt(fn, comments=['#'], usecols=(0,1), unpack=True)
ptd=0
if df > 0.0:
    ptd = int(len(x)*df)
if dp > 0:
    ptd=dp
x = x[ptd:]
y = y[ptd:]
popt, pcov = curve_fit(func, x, y)
av=np.mean(y)
sd=np.std(y)
plt.ylim(av-sd*sdf,av+sd*sdf)
plt.plot(x, y, '-', label='H(t)')
plt.plot(x, func(x, *popt), '-', label='fit: y=(%G)*x + (%G)' % tuple(popt))
plt.legend()
plt.show()
