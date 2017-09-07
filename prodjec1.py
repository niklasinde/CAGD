#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 14:17:05 2017

@author: niklasinde
"""

import numpy as np



class convexhull:
    def __init__(self,points,low = True):
        self.points = points
        
        self.upper = copy(points)
        self.lower = copy(points)
        
        


        
    
    def calclutning(self,p0,points,low):
        # Should we calculate the upper or lower points in the convex hull 
        # If self.low = True we calculate the upper points
        # p_0 calculates the angel for P_0 to the rest of the points
        
        a = []
        for i in range(len(points)):
            # If points[i] < p0 remove the points from the list
            if points[i] < p0:
                pass
            if points[i] == p0:
                
            if p[0][1]-p[i][1]== 0:
                b = (-np.inf)
            else:
                b = (p[0][0]-p[i][0])/(p[0][1]-p[i][1])
        a.append(b)
        return(a)
    

        
    def findpoints(self):
        pass
    
    
    
def DivDiff(p):



    
def convexhull(points):
    points.sort()
    
    upper = []
    lower = []
    lowerpoints = copy(points)
    upperpoints = copy(points)
    
    upper.append(lowerpoints[0])
    lower.append(upperpoints[0])
    del lowerpoints[0]
    del upperpoints[0]
    while True:
        lutning = DivDiff(points)
        maxindex = lutning.index(max(lutning))
#        print(points,lutning)
#        print(maxindex)
        upper.append(points[maxindex])
        deletebadx(points,points[maxindex])
        del points[maxindex]
        
        if len(points)==0:
            break
        
    return(upper,lower)
    
points= [[1,2],[31/6,3],[-8/3,4],[0,5]]
print(convexhull(points))
#print(points.sort())