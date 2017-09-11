#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 14:17:05 2017

@author: niklasinde
"""

import numpy as np
import matplotlib.pyplot as plt

class convexhull:
    def __init__(self,points):
        """ The points need to be in the form [[x_0,y_o],[x_1,y_1],[x_2,y_2]...[x_n,y_n]]"""
    
        self.points = points
        self.loader()
        
    
    def _upperpoints(self,p0,p):
        """ Finds the upper points of the convex hull"""
        yvalue = np.inf
        lutning = -np.inf
        results = None
        delete = []
        for i in range(len(p)):
            # if P_i,x < p_0,x
            if p[i][0] < p0[0]:
                delete.append(i)
            elif p[i][0]==p0[0]:
                if p0[1]< p[i][1] and p[i][1] < yvalue:
                    yvalue = p[i][1]
                    lutning = np.inf
                    results = p[i]
            else:
                dydx =(p0[1]-p[i][1])/(p0[0]-p[i][0])
                if lutning < dydx:
                    lutning = dydx 
                    results = p[i]
        for i in reversed(delete):
            del p[i]
#        del p[0]
        return (results,p)        
    
    def _lowerpoints(self,p0,p):
        """ Lower points we could put this function and upper points function together but this it for now. """
        yvalue = np.inf
        lutning = np.inf
        results = None
        delete = []
        for i in range(len(p)):
            
            # if P_i,x < p_0,x
            if p[i][0] < p0[0]:
                # We have to save the deleted index to a list otherwice the foreloop index gets messed up
                delete.append(i)  
            # if p_i,x == p_0,x
            elif p[i][0]==p0[0]:
                if  p0[1]> p[i][1] and p[i][1] < yvalue:
                    yvalue = p[i][1]
                    results = p[i]
                    lutning = -np.inf
            else:
                dydx =(p0[1]-p[i][1])/(p0[0]-p[i][0])
                if lutning > dydx:
                    lutning = dydx 
                    results = p[i]
        #delete x points that we dont need for next iteration.
        for i in reversed(delete):
            del p[i]
        return (results,p)        

        
    def _loader(self):
        """ function to calculate the upper and lower points of the cunvex hull """
        points = self.points
        points = sorted(points)
        
        upper = [points[0]]
        lower = [points[0]]
        
        pi , ptsu  = self._upperpoints(upper[-1],points[1:])
        pj, ptsl = self._lowerpoints(lower[-1],points[1:])
        while True:
            upper.append(pi)
            pi , ptsu = self._upperpoints(upper[-1],ptsu)
            if pi == None:
                break
        self.upper = upper
        
        while True:
            lower.append(pj)
            
            pj , ptsl = self._lowerpoints(lower[-1],ptsl)

            if pj == None:
                break
        # Connection the dots.
        if lower[-1][0]==upper[-1][0] and lower[-1][1]!=upper[-1][1]:
            lower.append(upper[-1])
        self.lower = lower
        
    def __call__(self,x:"list[x,y] or [[x_0,y_0],...,[x_n,y_n]])":       
        """A function to check if x is in the convex hull. 
            If x is in the convex hull.
            If x is a list it will return True 
        #If we put in a list of x,y values this is still possible 
        # Returns true if all of the the points are in the convexhull."""
        if  np.shape(np.array(x))!=(2,):
            X = [self.checker(x[i]) for i in range(len(x))]
            if X==True:
                return(True)
            else:
                return(False)
        low = self.lower
        up = self.upper
        
        
        if not self.upper[0][0] <= x[0] <= self.upper[-1][0]:
            return False
        
        
        if x[0] == up[0][0] and up[0][0] == up[1][0]:
            if up[0][1]<= x[1] <= up[1][1]:
                return(True)
            else:
                return(False)
                
                
        elif  x[0] == low[-1][0] and low[-1][0]==[-2][0]:
            if low[-2][1]<= x[1] <= low[-1][1]:
                return(True)
            else:
                return(False)
        else:
            for i in range(len(up)-1):
                if up[0][0] <= x[0] <= up[i+1][0]:
                    up_index = i 
                    
            for i in range(len(low)-1):

                if low[i][0] <= x[0] <= low[i+1][0]:

                    low_index = i 
            
            # calculate the function of the line between the two dots which are found in the function above.
            klower =(low[low_index][1]-low[low_index+1][1])/(low[low_index][0]-low[low_index+1][0])
            kupper =(up[up_index][1]-up[up_index+1][1])/(up[up_index][0]-up[up_index+1][0])
            
            mlower = low[low_index][1]-klower*low[low_index][0]
            mupper = up[up_index][1]-kupper*up[up_index][0]
            if klower*x[0]+mlower < x[1] < kupper*x[0]+mupper:
                return(True)
            else:
                return(False)
        
    def plot(self):
        low = self.lower
        up = self.upper
        points = self.points
        x = [points[i][0] for i in range(len(points))]
        y = [points[i][1] for i in range(len(points))]
        plt.plot(x,y ,"o")
        
        x1,x2,y1,y2 = up[0][0],up[1][0],up[0][1],up[1][1]
        upline, = plt.plot((x1,x2),(y1,y2),'-',label="Upper Limit",color="red")
        for i in range(1,len(up)-1):
            x1,x2,y1,y2 = up[i][0],up[i+1][0],up[i][1],up[i+1][1]
            upline, = plt.plot((x1,x2),(y1,y2),'-',color="red")

        x1,x2,y1,y2 = low[0][0],low[1][0],low[0][1],low[1][1]
        lowline, = plt.plot((x1,x2),(y1,y2),'-',label="Lower Limit",color="green")
        for i in range(1,len(low)-1):
            x1,x2,y1,y2 = low[i][0],low[i+1][0],low[i][1],low[i+1][1]
            lowline, = plt.plot((x1,x2),(y1,y2),'-',color="green")
        plt.legend()
        plt.title("Convex hull of given points.")
        plt.show()


    
points = [[0.05, 0.02], [0.1, 0.2],[1,1]]

#x = [points[i][0] for i in range(len(points))]
#y = [points[i][1] for i in range(len(points))]
#plt.plot(x,y ,"o")
a = convexhull(points)
a.plot()
print(a.checker([[0.7,0.2],[0.9,-0.1],[1.3,0]]))
print(a.checker([0,0]))


