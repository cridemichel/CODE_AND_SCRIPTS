#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 16:44:21 2019

@author: demichel
"""
import numpy as np
class mat:
    def __init__(self):
        self.A=[[1,2],[3,4]]
    A = []
    def __add__(self,B):
        D=[]
        for i in range(0,len(self.A)):
            C=[]
            for j in range(0,len(self.A[i])):
                 C.append(self.A[i][j]+B.A[i][j])
            D.append(C)
        return D
#A=[[1,2,3],[4,10,6],[7,8,9]]
B=[[1,2,3],[4,5,6],[7,8,9]]
F=B[1:3][1][-1]
A=mat()
B=mat()
print ('A=',A.A, ' len(A)=',len(A.A))
print ('A+B=',A+B)
quit()
print('boh=',B[1:3])
print ('>>>F=',F)
quit()
M=B
A=('a','b','c')
c=M[0]
print (c)
d=[M[i] for i in range(0,3)]
print(d)
a=1,2,3
C=[row[1] for row in M if row[1] % 2 == 0] 
print (C)
x=3.2
print(np.arccos(0.321))
#import itertools
#from itertools import izip
#print(A+B)
print("tupla=",(2,)+A[1:len(A)])
D={'a': 'mario', 'd': 'pino', 'b':'luigi'}
LD=list(D.keys())
LD.sort()
print("sorted=",LD)
if 'a' in D:
    print(D['a'])
T1={'a','b','c'}
T2={'b','c'}
T3=T1 | T2
print('T3prima=',T3)
L1=list(T3)
print('zip(L1)=',zip(T1))
print('L1=',L1)
L1.sort()
print('L1sorted=',set(L1))
T3=set(L1)
print('T3=',T3)
L1=[1,4,53,2,-10,120,2,4]
print("mas(L1)=",max(L1))