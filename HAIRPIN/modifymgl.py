#!/usr/local/bin/python3
import sys
args=sys.argv
with open(args[1]) as f:
    lines=f.readlines()
L=[ 55.71817679104672294, 78.050693589422294849, 102.33563723371483434 ]
L2=L[2]/2.0
L3=L[2]/3.0
newf=[]
cc=0
colors=['DeepSkyBlue', 'turquoise', 'orange','red']
for l in lines:
    p=l.strip('\n').split()
    layer=int((float(p[2])+L2)/L3)
    if (cc % 2 == 0):
        p[9]='C['+colors[layer]+']'
        newf.append(p)
    else:
        p[4]='0.5'
        p[5]='C['+colors[3]+']'
    s=' '.join(p)
    print(s)
    cc+=1
