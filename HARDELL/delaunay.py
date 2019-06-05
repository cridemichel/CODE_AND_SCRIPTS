#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 09:28:48 2019

@author: demichel
"""
import sys
#import numpy as np
from numpy import loadtxt
import matplotlib.pyplot as plt
#from scipy.spatial import Voronoi, voronoi_plot_2d
#import numpy as np
from scipy.spatial import Delaunay
#restituisce un array con tutti i primi vicini del vertice pindex (vedi sotto)
def find_neighbors(pindex, triang): 
    return triang.vertex_neighbor_vertices[1][triang.vertex_neighbor_vertices[0][pindex]:triang.vertex_neighbor_vertices[0][pindex+1]] 
points = loadtxt("/Users/demichel/conf-python", comments="#", delimiter=" ", unpack=False)
#points = np.array([[0, 0], [0, 1.1], [1, 0], [1, 1]])
#vor=Voronoi(points)
tri = Delaunay(points)
#print (points)
#print ("number of points=",points.size/2)
#for pindex in range(0,2515):
#    nn=find_neighbors(pindex,tri)
#    print("i=",pindex, nn)
#plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
#plt.plot(points[:,0], points[:,1], 'o')
#plt.show()
#sys.exit(0)
#voronoi_plot_2d(vor)
#plt.plot(points[:,0], points[:,1], 'o')
#plt.show()
#nnmy = loadtxt("/Users/demichel/nnmy.dat", dtype="int",comments="#", delimiter=" ", unpack=False)
filename="/Users/demichel/nnmy.dat"
#nnmy=np.genfromtxt(filename, delimiter=" ")
#print(filecontents)
pindex=0
with open(filename,"r") as f:
    for line in f:
        nnmy=line.strip('\n').split(' ')    
        print ("checkin i=", pindex)
        nn=find_neighbors(pindex,tri)
        if len(nn) != len(nnmy):
            print ("WARNING numero diverso elementi per pindex=", pindex)
            print ("nn=", nn)
            print ("nnmy=", nnmy)
            sys.exit(0)
            print("len=", len(nn))
            for p2 in range(0,len(nn)):
                trovato=0         
                for p3 in range(0,len(nnmy)):
                    print ("nn=", nn[p2], " nnmy=", nnmy[p3])
                    if int(nn[p2]) == int(nnmy[p3]):
                        print ("qui\n")
                        trovato=1
                        break
                    if trovato == 0:
                        print("WARNING non trovato\n")
                        sys.exit(0)
        pindex += 1
print ("tutto ok!\n")                