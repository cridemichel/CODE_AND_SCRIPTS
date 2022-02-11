#!/usr/local/bin/python3
import numpy as np, sys
if len(sys.argv) < 4:
    print('You must supply name of the file and min and max abscissa')
fn=sys.argv[1]
xmin=float(sys.argv[2])
xmax=float(sys.argv[3])
[x, y] = np.loadtxt(fn, usecols=(0,1), delimiter=' ', unpack=True)
sy=[] # type: list[float]
cc=0
for el in x:
    if el <= xmax and el >= xmin:
        sy = np.append(sy, y[cc])
    cc = cc+1
va = np.var(sy)
m = np.mean(sy)
n = len(sy)
print('Average: ', m, ' StdDev=', va,  ' #el=', n)
