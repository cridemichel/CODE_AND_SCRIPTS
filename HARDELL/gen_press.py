#!/usr/bin/env python3
import os,sys
start=float(sys.argv[1])
step=float(sys.argv[2])
n=int(sys.argv[3])
for i in range(n):
    print(start+step*i,' ', end='')
print('')
