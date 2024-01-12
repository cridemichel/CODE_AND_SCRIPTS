import sys
import numpy as np
args=sys.argv
with open(args[1], "r", encoding='utf8') as f1:
    cA=f1.readlines()
args=sys.argv
with open(args[2], "r", encoding='utf8') as f2:
    cB=f2.readlines()
cc=0
somma=0.0
sc=0.0
for l in cA:
    if cc >= 1:
        #print("l=", l)
        l1 = l.split()
        l2 = cB[cc].split()
        if len(l1) == 3 and len(l2)==3:
            rA=np.array([float(l1[0]),float(l1[1]), float(l1[2])])
            rB=np.array([float(l2[0]),float(l2[1]), float(l2[2])])
            rAB=rA-rB
            somma+=np.dot(rAB,rAB)
            sc=sc+1.0
    cc=cc+1
somma = somma/float(sc)
print('difference is:', np.sqrt(somma))

