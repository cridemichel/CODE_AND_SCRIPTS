#!/usr/local/bin/python3
import numpy as np, sys
import matplotlib.pyplot as plt
if len(sys.argv) == 1:
    print('You must supply a file with data')
fn=sys.argv[1]
if len(sys.argv) > 2:
    maxnb=int(sys.argv[2])
else:
    maxnb=1000
x, y = np.loadtxt(fn, usecols=(0,1), delimiter=' ', unpack=True)
va = np.var(y)
m = np.mean(y)
n = len(y)
xba=[] # type: list[float]
yba=[] # type: list[float]
for t in range(2, maxnb):
    nblocks = int(n/t)
    sbl=np.array_split(y, nblocks)
    #if t == maxnb-1:
    #    print('sbl=', sbl)
    vb=0.0
    cc=0
    blav=[] # type: list[float]
    for v in sbl:
        cc = cc + 1
        blav = np.append(blav, np.var(v))
    vb = np.var(blav)
    #if int(t / 100) == t/100:
    print('nblocks=', nblocks, 'len(sbl)=', len(sbl), 't=', t, ' vb=', vb, 'v=', va)
    xba=np.append(xba,t)
    yba=np.append(yba,t*vb/va)
np.savetxt('blockavg.dat', np.column_stack([xba,yba]), fmt='%.18g', delimiter=' ')
plt.plot(xba, yba, 'ob')
print('Error=', np.sqrt(va*yba[-1]/n))
plt.show()
