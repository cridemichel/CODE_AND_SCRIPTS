#!/usr/bin/env python3
import sys, os
arg=sys.argv
if len(arg) < 4:
    print('syntax is:\n set_one_par.py <parfile> <parname> <value>')
    quit()
par=arg[2]
val=arg[3]
with open(arg[1]) as f:
	lines=f.readlines()
nl=[]
for l in lines:
	l0=l.strip('\n').split('#')
	comment='#'+'#'.join(l0[1:])
	la=l0[0].split(':')
	if la[0].strip(' ')==par:
		nl.append(la[0]+': '+val+' '+comment+'\n')
	else:
		nl.append(l)	

with open('_aaa_', 'w') as f:
	for l in nl:
		f.write(l)
os.system('mv _aaa_ '+arg[1])
