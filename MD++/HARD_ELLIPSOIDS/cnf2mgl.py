#!/usr/bin/python3
import os,sys
if len(sys.argv)==1:
    print('please supply a filename')
    quit
with open(sys.argv[1]) as f:
    lines=f.readlines()
cc=0
for l in lines:
    if cc > 0:
        print(l.strip('\n')+' @ 1.4 0.5 0.5 C[red]')
    c++
