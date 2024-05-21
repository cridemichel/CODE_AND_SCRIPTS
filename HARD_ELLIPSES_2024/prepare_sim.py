#!/usr/bin/python3
import os, sys
import numpy as np
curd=os.getcwd()
d=os.path.split(curd.strip('\n'))[1]
args=sys.argv
#if len(args)>1:
#    T=float(args[1])
#else:
#    T=0.125
q=float(d.split('_')[1].strip(' '))
print ('q=',q)
os.system('cp ../phscpars_templ.xasc .')
lines2new=[]
with open("phscpars_templ.xasc") as f:
    lines2=f.readlines()
    for l2 in lines2:
        #print ('l2=', l2)
        ls=l2.strip('\n').split(' ')
        ln = ' '.join([l.replace('_Q_',str(q)) for l in ls])+'\n'
        #ln = ln.replace('_A_',sa)
        lines2new.append(ln)
with open("phscpars.xasc","w") as f:
    for l in lines2new:
        #print ('l=',l)
        f.write(l)
