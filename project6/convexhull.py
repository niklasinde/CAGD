#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 14:17:05 2017

@author: niklasinde
"""

import numpy as np
import matplotlib.pyplot as plt

class convexhull:
    def __init__(self,points:"points"):
        """ The points need to be in the form [[x_0,y_o],[x_1,y_1],[x_2,y_2]...[x_n,y_n]]"""
        self.points = points
        self._loader()
        self.domain = [self.lower[0][0],self.lower[-1][0]]
        
    
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
        if type(self.points)== np.ndarray:
#            print("yes")
            self.points = self.points.tolist()
        points = sorted(self.points)
        
        
        upper,lower = [points[0]],[points[0]]
        ptsu,ptsl  = points[1:],points[1:]
        
        while True:
            pi , ptsu = self._upperpoints(upper[-1],ptsu)
            if pi == None:
                break
            else:
                upper.append(pi)
        self.upper = upper
        
        
        while True:
            pj , ptsl = self._lowerpoints(lower[-1],ptsl)

            if pj == None:
                break
            else:
                lower.append(pj)
                
                
        # Connection the dots.
        if lower[-1][0]==upper[-1][0] and lower[-1][1]!=upper[-1][1]:
            lower.append(upper[-1])
        self.lower = lower   
    def __call__(self,x:"list[x,y] or [[x_0,y_0],...,[x_n,y_n]])",ret = False):       
        """A function to check if x is in the convex hull. 
            If x is in the convex hull.
            If x is a list it will return True 
        #If we put in a list of x,y values this is still possible 
        # Returns true if all of the the points are in the convexhull."""
        if  np.shape(np.array(x))!=(2,):
            X = [self.__call__(x[i]) for i in range(len(x))]
            if ret:
                return X
            else:
                if X==True:
                    return(True)
                else:
                    return(False)
                        
        low = self.lower
        up = self.upper
        
        if not self.lower[0][0] <= x[0] <= self.lower[-1][0]:
#            print(1)
            return False
        
        
        if x[0] == up[0][0] and up[0][0] == up[1][0]:
            if up[0][1]<= x[1] <= up[1][1]:
                
                return(True)
            else:
                print(2)
                return(False)
                
                
        elif  x[0] == low[-1][0] and low[-1][0]==low[-2][0]:
            if low[-2][1]<= x[1] <= low[-1][1]:
                return(True)
            else:
#                print(2)
                return(False)
        else:
            for i in range(len(up)-1):
                if up[0][0] <= x[0] <= up[i+1][0]:
                    up_index = i 
                    
            for i in range(len(low)-1):

                if low[i][0] <= x[0] <= low[i+1][0]:

                    low_index = i 
            if low[low_index][0]==low[low_index+1][0]:
                raise ZeroDivisionError(str(low),str(x))
            if up[up_index][0]==up[up_index+1][0]:
                raise ZeroDivisionError(str(up)+str(x))
            # calculate the function of the line between the two dots which are found in the function above.
            klower =(low[low_index][1]-low[low_index+1][1])/(low[low_index][0]-low[low_index+1][0])
            kupper =(up[up_index][1]-up[up_index+1][1])/(up[up_index][0]-up[up_index+1][0])
            
            mlower = low[low_index][1]-klower*low[low_index][0]
            mupper = up[up_index][1]-kupper*up[up_index][0]
            if klower*x[0]+mlower <= x[1] <= kupper*x[0]+mupper:
                return(True)
            else:
#                print(3)
                return(False)
    def area(self):
        """
        Calculating the area of the points using shoelace formula.
        https://en.wikipedia.org/wiki/Shoelace_formula
        """
        x = [k[0] for k in self.points]
        y = [k[1] for k in self.points]
        return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
    def intersect(self,other):
        # Other == line 
        if other.points[0][0]==other.points[1][0]: raise ValueError("Line is vertical, this is coming soon." )
        else:
            k =(other.points[0][1]-other.points[1][1])/(other.points[0][0]-other.points[1][0])
            m = other.points[0][1]-k*other.points[0][0]
            f = lambda x: k * x + m
            for i in range(len(self.points)-1):
                if np.sign(f(self.points[i][0])-self.points[i][1]) != np.sign(f(self.points[i+1][0])-self.points[i+1][1]):
                    return(True)
            return False
        
    def render(self):
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
    
    def plot(self):
        self.render()
        plt.legend()
        plt.title("Convex hull of given points.")
        plt.show()


    
#points = [[0.05, 0.02], [0.1, 0.2],[1,1]]

#x = [points[i][0] for i in range(len(points))]
#y = [points[i][1] for i in range(len(points))]
#plt.plot(x,y ,"o")
#a = convexhull(points)
#a.plot()
#print(a([[0.7,0.2],[0.9,-0.1],[1.3,0]]))
#print(a([0,0]))


