#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 09:28:48 2019

@author: demichel
"""
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
print (points)
print ("number of points=",points.size/2)
for pindex in range(0,2515):
    nn=find_neighbors(pindex,tri)
    print("i=",pindex, nn)
plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
plt.plot(points[:,0], points[:,1], 'o')
plt.show()
#voronoi_plot_2d(vor)
plt.plot(points[:,0], points[:,1], 'o')
plt.show()