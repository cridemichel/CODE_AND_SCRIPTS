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
fn=argv[1]
if len(argv) > 2:
    sdf = float(argv[2])
else:
    sdf = 30.0
f, axarr = plt.subplots()
axarr.tick_params(labeltop=True, labelright=True)
plt.title(r'orthoterphenyl $\rho\approx 1.106 g/cm^3$ $ T=281 K $ N=2352 $\delta t=5 fs$')
plt.xlabel('t (ps)')
plt.ylabel('energy')
x, y = np.loadtxt(fn, comments=['#'], usecols=(0,1), unpack=True)
popt, pcov = curve_fit(func, x, y)
av=np.mean(y)
sd=np.std(y)
plt.ylim(av-sd*30,av+sd*sdf)
plt.plot(x, y, '-', label='H(t)')
plt.plot(x, func(x, *popt), '-', label='fit: y=(%G)*x + (%G)' % tuple(popt))
plt.legend()
plt.show()
