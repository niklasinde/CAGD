#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 19:03:11 2017
 
@author: niklasinde
"""
import sys
 


import numpy as np
from numpy import *
from sympy.utilities.lambdify import lambdify
import sympy as sm
import matplotlib.pyplot as plt
import sympy
import scipy
plt.rcParams["figure.figsize"] = [10,10]
import time as time
 
class spline():
    """ Class to calculate a spline"""
    def __init__(self, knot:"knotvector", k = 3, coeff = None):
        if len(knot)<= k+1:
            raise SyntaxError("We need atleast k+1 number of points in the knotvector")
        self.k = k
        self.points = knot
        if not type(coeff) == np.ndarray:
            self.coeff = [1 for i in range(len(self.points)-self.k-1)]
        else:
            self.coeff = coeff
        self.x = sm.Symbol("x")
#        self.lastx = kn
        self.N,self.F =self.basicfunction()
 
    #### Loading functions ###
    """ functions that will load when a variable is called with this class"""
    def ref(self, p):
        """creates a list with values between the x values in points
        because sympy.Symbol("x") cant select a value and we need a x-value
        in the recurrence function"""
        return [p[i]+(p[i+1]-p[i])/2 for i in range(len(p)-1)]
     
    def getrelpts(self):
        n = len(self.points)-self.k-1
        p= [self.points[i:i+self.k+2] for i in range(n)]
        return(p)
    def movepoint(self,new_point,index_oldpoint):
        updatebasislist = [index_oldpoint]
           
    def updatebasis(self, newpoints):
        N1,F1 = self.N,self.F
        self.points = newpoints
        N2,F2 = self.basicfunction()
        
        
    def rec(self, p, i, k, xi):
        """recurrence formula for b-spline"""
        # We should lambdify this function for higer speed.
        if k == 1:
            if p[i] == p[i+1]:
                return 0
            if  p[i] <= xi <= p[i+1]:
                return 1
            else:
                return 0
        else:
            def div(lhs, rhs):
                if rhs == 0:
                    return 0
                else:
                    return lhs / rhs
            u =  (div((self.x-p[i]),(p[i+k-1]-p[i]))*self.rec(p,i,k-1,xi)+
                  div((p[i+k]-self.x),(p[i+k]-p[i+1]))*self.rec(p,i+1,k-1,xi))
            return u

            
        
    def basicfunction(self):
        n1 = []
        f1 = []
        p = self.getrelpts()
         
        for i in range(len(p)):
            n2=[]
            f2 = []
            xi = self.ref(p[i])
            for j in range(self.k+1):
                func = self.rec(p[i],0,self.k+1,xi[j])
                evalfunc = lambdify(self.x, func, modules=['numpy'])
                n2.append(evalfunc)
                f2.append(func)
            n1.append(n2)
            f1.append(f2)
             
        return n1,f1
 
    def basisplot(self):
        """ Plots all the basis functions """
        for i in range(len(self.N)):
            for j in range(len(self.N[i])):
                if self.N[i][j]!= 0:
                    x = linspace(self.points[i+j],self.points[i+j+1],50)
                    y = self.coeff[i]*self.N[i][j](x)
                    plt.plot(x,y)
        plt.title("Plot of the basic functions for the splines")
#        plt.show()
 
    def evalfull(self,X):
        p = self.points
        if np.shape(np.array(X)) != ():
            self.test = True
            return [self.evalfull(x) for x in X]
         
        hot_interval = self.hotinterval(X)
        func_val = 0
        for i in range(len(self.N)):
            if hot_interval < i or i+self.k < hot_interval:
                pass
            else:
                evalfunc = self.N[i][hot_interval-i](X)
                func_val += self.coeff[i]*evalfunc
        return func_val
     
     
    def evalvector(self,control,X):
        p = self.points
        Y = []
        y = array([0,0])
#        X =np.linspace(self.p[0],self.p[-1])
#        if len
        for x in X:
            hi = self.hotinterval(x)
            func_val = np.array([0.,0.])
            for i in range(len(self.N)):
                if hi < i or i+ self.k < hi:
                    pass
                else:
                    evalfunc = self.N[i][hi-i](x)
                    func_val += control[i]*np.array(evalfunc)
            Y.append(func_val)
        Y = np.array(Y)
        return Y
     
             
    def evalbasis(self,x,i):
        p = self.points
        func_val=0
         
        hot_interval = self.hotinterval(x)
        if hot_interval <= i or i+self.k <= hot_interval:
            return 0.
        else:
            evalfunc = self.coeff[i]*self.N[i][hot_interval-i](x)
            return evalfunc
    def evalbasisi(self,i,x):
        return self.N[i](x)
     
    def hotinterval(self, x):
        p = self.points
        for i in range(len(p)-1):
            if p[i] <= x and x <= p[i+1]:
                return i
        raise ValueError(x," is not in the interval")
         
 
    def hotinterval_(self, t):
        if not self.points[0] <= t <= self.points[-1]:
            raise ValueError("t is not in the interval.")
             
        pts = self.points
        first,mid,last = 0, last // 2, len(pts) - 1
        u_left, u_right = u[mid],  u[mid + 1]
         
        while (not (u_left <= t <= u_right)): # Note: <= on both sides, since nodes might be coincident
            if (t > u_left):
                first = mid + 1
            else:
                last = mid
            mid = (last + first) // 2
            u_left = u[mid]
            u_right = u[mid + 1]
        return mid
     
class matrixequation:
    """ Calculates A d = x  """
 
    def __init__(self,xy,points):
        self.xy = xy
        self.p = points
        self.N = spline(points)
 
        self.xyvec = self.xyvector()
        self.e = self.elist()
        self.A = self.Avector()

        self.z = self.solver()
        
    def Avector(self):
        """creats a square tridiagonal matrix of size n"""
        n = len(self.N.N)
        A = np.zeros((n,n))
        for col in range(n):
            for row in range(n):
                A[row, col] = self.N.evalbasis(self.e[row],col)
        return A
 
    def elist(self):
        l = []
        p = self.p
        for i in range(1,len(self.N.N)+1):
            l.append((p[i]+p[i+1]+p[i+2])/3)
        return l
 
    def xyvector(self):
        """ When i write xy its because its the same function for x and y"""
        xy = np.array(self.xy)
        xy.reshape((1,len(self.xy)))
        return xy
 
    def solver(self):
        print(np.shape(self.xy),np.shape(self.A))
        print("HWJ")
        z = scipy.linalg.solve(self.A,np.array(self.xy).transpose())
        self.coeff = z
        return(z)
 
 
class interpolation:
    def __init__(self,points,x,y, k = 3):
        self.k = k
        self.interx = x
        self.intery = y
        self.pts = points
        self.xcoeff = matrixequation(x,self.pts)
        self.ycoeff = matrixequation(y,self.pts)
        
        self.splinex = spline(self.pts,self.k,self.xcoeff.coeff)
        
        self.spliney = spline(self.pts,self.k,self.ycoeff.coeff)
    def x(self,x):
        return(self.splinex.evalfull(x))
    def y(self,y):
        return(self.spliney.evalfull(y))
    def plotbasis(self):
        print("xbasis:")
        self.splinex.basisplot()
        print("ybasis:")
        self.spliney.basisplot()
    def plotinter(self,pts = 1000):
        t0 = time.time()
        dt = 1/pts
        t = arange(self.pts[0], self.pts[-1] + dt, dt)
        plt.plot(self.x(t),self.y(t))
        plt.scatter(self.interx,self.intery,color="red")
        print(time.time()-t0)
 

knot = [0,0,0,0,1,2,3,4,5,6,7,7,7,7]
knot2 = [0,0,0,0,1,2,3.9,4,5,6,7,7,7,7]

s = spline(knot)
s.basisplot()
#plt.show()
s1 = spline(knot2)
s1.basisplot()



