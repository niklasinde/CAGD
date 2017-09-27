#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 19:03:11 2017

@author: niklasinde
"""
import sys

#plt.rcParams["figure.figsize"] = [6,3]

import numpy as np
from numpy import *
import matplotlib as plt
from sympy.utilities.lambdify import lambdify
import sympy as sm
import matplotlib.pyplot as plt
import sympy
import scipy

import time as time

class spline():
    """ Class to calculate a spline"""
    def __init__(self, points:"knotvector", k = 3, coeff = None):
        print(len(points),k+1)
        if len(points)<= k+1:
            raise SyntaxError("We need atleast k+1 number of points in the knotvector")
        self.k = k
        self.points = points
        if not type(coeff) == np.ndarray:
            self.coeff = [1 for i in range(len(points)-self.k-1)]
        else:
            self.coeff = coeff
#        self.k = deg
        self.x = sm.Symbol("x")
        self.loader()
        self.lastx = points
    #### Loading functions ###
    """ functions that will load when a variable is called with this class"""
    
    def ref(self, p):
        """creates a list with values between the x values in points
        because sympy.Symbol("x") cant select a value and we need a x-value
        in the recurrence function"""
        xi =[p[i]+(p[i+1]-p[i])/2 for i in range(len(p)-1)]
        return xi
    def getrelpts(self):
        n = len(self.points)-self.k-1
        p = [self.points[i:i+self.k+2] for i in range(n)]
        return(p)
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
#                print(func)
            n1.append(n2)
            f1.append(f2)
        self.N = n1
        self.F = f1

    def loader(self):
        ### loading function ###
        ### Call this if the points are updated ###
        self.basicfunction()

    def basisplot(self):
        """ Plots all the basis functions """
        for i in range(len(self.N)):
            for j in range(len(self.N[i])):
                if self.N[i][j]!= 0:
                    x = linspace(self.points[i+j],self.points[i+j+1],50)
                    y = self.coeff[i]*self.N[i][j](x)
                    plt.plot(x,y)
        plt.title("Plot of the basic functions for the splines")
        plt.show()

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
    
    
    def evalvector(self,vector, X):
        p = self.points
        Y = []
        y = array([0,0])
        for x in X:
            hi = self.hotinterval(x)
            func_val = np.array([0.,0.])
            for i in range(len(self.N)):
                if hi < i or i+ self.k < hi:
                    pass
                else:
                    evalfunc = self.N[i][hi-i](x)
                    func_val += vector[i]*np.array(evalfunc)
#                    print(func_val,np.array(evalfunc),"\n",vector[i])

            Y.append(func_val)
        Y = np.array(Y)
        return Y
    
            
    def evalbasis(self,x,i):
        p = self.points
        func_val=0
        hot_interval = self.hotinterval(x)
        # If the basis is 0 at x
        if hot_interval <= i or i+self.k <= hot_interval:
            return 0.
#        if xi is in the function value return 0
        else:
            evalfunc = self.coeff[i]*self.N[i][hot_interval-i](x)
            return evalfunc
    def hotinterval(self, x):
        p = self.points
        for i in range(len(p)-1):
            if p[i] <= x and x <= p[i+1]:
                return i
        raise ValueError(x," is not in the interval")
        
    def evalbasisi(self,i,x):
        return self.N[i](x)
    def hotinterval_(self, t):
        if not self.points[0] <= t <= self.points[-1]:
            raise ValueError("t is not in the interval.")
        """
        Binary search to find hotspot for t.
        
        input:
            t: one value somewhere in the knot vector.
        
        Output: 
            Indicies for the knot vector which are precisely smaller and larger in the padded knot vector.
        """

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
class matrixequation():
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
        z = scipy.linalg.solve(self.A,np.array(self.xy).transpose())
        self.coeff = z
        return(z)


class interpolation():
    def __init__(self,points,x,y, k = 3):
        self.k = k
        self.interx = x
        self.intery = y
        self.pts = points
        self.xcoeff = matrixequation(x,self.pts)
        self.ycoeff = matrixequation(y,self.pts)
        print(self.pts,"pts")
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


class approximation():
    def __init__(self,pts,xy):
        if type(xy) == list:
            self.xy = np.array(xy)
        else:
            self.xy = xy
        self.pts = pts
        
        spl = spline(pts,3,xy)
        
        
        
        
### Task 1 
        
#knot = [0,0,1,1]
#A = spline(knot,k=1)
#A.basisplot()
#print(A.F)
#knot2 = [0,1,2,3]
#A = spline(knot2,k=1)
#A.basisplot()
#print(A.F)
#knot2 = [0,0,0,1,1,1]
#A = spline(knot2,k=2)
#A.basisplot()
#print(A.F)
#knot2 = [0,1,2,3,4,5]
#A = spline(knot2,k=2)
#A.basisplot()
#print(A.F)

#x = linspace(0,1,50)
#plt.plot(x,A.evalfull(x))
### task 2


#knot3 =[0,0,0,0.3,0.5,0.5,0.6,1,1,1]
#A = spline(knot3, k= 3)
#A.basisplot()

def hotinterval(knot,x):
    p = knot
    for i in range(len(p)-1):
        if p[i] <= x and x <= p[i+1]:
            return i
    raise ValueError(x," is not in the interval")
        
def hotInterval(knot,t):
    if not knot[0] <= t <= knot[-1]:
        raise ValueError("{} is not in the interval".format(t))
    first,last = 0, len(knot) - 1
    mid = last // 2
    k_left, k_right = knot[mid],  knot[mid + 1]
    while (not (k_left <= t <= k_right)):
        if (k_left<t):
            first = mid + 1
        else:
            last = mid
        mid = (last + first) // 2
        k_left = knot[mid]
        k_right = knot[mid + 1]
    if knot[mid] == knot[-1]:
        mid = mid-1
    return mid
def fullsupport(knot, t, deg):
    return knot[deg] <= t <= knot[-deg-1]

knot = [0, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1]
pts = [0, 0.12, 0.24, 0.4, 0.53, 0.78, 0.8, 1]
for i in pts:
    print(i,fullsupport(knot,i,2))
a = spline(knot,k=2)

x= linspace(knot[0],knot[-1],100)
y = a.evalfull(x)
plt.plot(x,y)
plt.plot((0.1,0.1),(0,1))
plt.plot((0.9,0.9),(0,1))
plt.ylim(-0.1,1.1)
plt.plot()
plt.savefig("task4.pdf")
## Task 4
#
#knot = [0, 0, 0, 0.3, 0.5, 0.5, 0.6, 1, 1, 1]
#control = np.array([[0,0],[3,4],[7,5],[9,2],[13,1],[10,1],[7,1]])
#
#A = spline(knot,2)
#
#for i in range(len(knot)-1):
#    X = linspace(knot[i],knot[i+1],20)
#    Y = A.evalvector(control,X)
#    plt.plot(Y[:,0],Y[:,1])
#    plt.scatter(control[:,0],control[:,1])
#plt.scatter(6,3,s=300)
#plt.scatter(6.1,3,s=50,color="black")
#plt.plot((7,6),(1,0))
#plt.grid()    
#plt.savefig("interpol2.pdf")

#plt.show()
#
#x = [x[0] for x in control]
#y = [x[1] for x in control]
#interpolation(knot,x,y,k=3)


#interpolation()


