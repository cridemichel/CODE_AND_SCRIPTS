#!/usr/bin/env python3
import os, sys
fp='lista_press'
if os.path.exists(fp):
    with open(fp) as f:
        l=f.readline()
    P=l.strip('\n').split(' ')
else:
    P=[ '1.0', '1.5', '2.0', '2.5', '3.0', '4.0', '4.5','5.0', '5.5', '6.0', '6.5', '7.0', '7.5', '8.0', '8.5', '9.0', '9.5','10.0' ]
for p in P:
    if p!= '' and (not os.path.exists('P_'+p)):
        os.system('mkdir P_'+p)
