#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 19:03:11 2017

@author: niklasinde
"""
import sys

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
#        if len(points) > deg-2:
#            raise IndexError("You dont have enough points to mach the desired polynomial degree.\n Add points or lower the degree.")
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
        self.test = False
    #### Loading functions ###
    """ functions that will load when a variable is called with this class"""

    def ref(self, p):
        """creates a list with values between the x values in points
        because sympy.Symbol("x") cant select a value and we need a x-value
        in the recurrence function"""
#        p = self.points
        xi =[p[i]+(p[i+1]-p[i])/2 for i in range(len(p)-1)]

        self.e =[(p[i]+p[i+1]+p[i+2])/2 for i in range(len(p)-self.k)]
        self.xi = xi
#        return xi



    def rec(self, p, i, k, xi):
        """recurrence formula for b-spline"""
        # We should lambdify this function for higer speed.
        
        if k == 1:
            if p[i] == p[i+1]:
                return 0

            if  p[i] <= xi and xi <= p[i+1]:
                return 1

            else:
                return 0
        else:
            def div(lhs, rhs):
                if rhs == 0:
                    return 0
                else:
                    return lhs / rhs
            u =  div((self.x-p[i]),(p[i+k-1]-p[i]))*self.rec(p,i,k-1,xi)+div((p[i+k]-self.x),(p[i+k]-p[i+1]))*self.rec(p,i+1,k-1,xi)
            return u


    def basicfunction(self):
        n1 = []
        f1 = []
        for i in range(len(self.points)-self.k-1):
            n2=[]
            p = self.points[i:i+self.k+2]
            xi = self.ref(p)
            for j in range(self.k+1):
                func = self.rec(p,0,self.k+1,xi[j])
                evalfunc = lambdify(self.x, func, modules=['numpy'])
                n2.append(evalfunc)
            n1.append(n2)
        self.N = n1
        
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
            if self.test == True:
                raise ValueError("Fel i evalful")
            self.test = True

            return [self.evalfull(i) for i in X]
        hot_interval = self.hotinterval(X)
        func_val = 0 
        for i in range(len(self.N)):
            if hot_interval < i or i+self.k < hot_interval:
                pass
            else:
                evalfunc = self.N[i][hot_interval-i](X)
                func_val += self.coeff[i]*evalfunc
        self.test = False
        return func_val

    def evalbasis(self,x,i):
        p = self.points
        func_val=0
        hot_interval = self.hotinterval(x)
        # If the basis is 0 at x
        if hot_interval < i or i+self.k < hot_interval:
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
    def hotinterval_(self, t):
        """
        Binary search to find hotspot for t.
        
        input:
            t: one value somewhere in the knot vector.
        
        Output: 
            Indicies for the knot vector which are precisely smaller and larger in the padded knot vector.
        """
        
        def shift(mid): # shift index because of index.
            return (mid + self.k, mid + self.k + 1)
    
        u = self.points
        first = 0
        last = len(u) - 1
        mid = last // 2
        if u[mid] == u[-1]:
                return shift(mid)
        
        u_left = u[mid]
        u_right = u[mid + 1]
        
        while (not (u_left <= t <= u_right)): # Note: <= on both sides, since nodes might be coincident
        
            if (t > u_left):
                first = mid + 1
            else:
                last = mid
                
            mid = (last + first) // 2
            
            if u[mid] == u[-1]: # The last point corresponding to Ui+1.start is missing (t[-1]), and we must ensure that u_right doesnt go out of bounds
                return shift(mid)
            
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
    def __init__(self,points,x,y):
        self.interx = x
        self.intery = y
        self.pts = points
        self.xcoeff = matrixequation(x,self.pts)
        self.ycoeff = matrixequation(y,self.pts)
        self.splinex = spline(points,self.xcoeff.coeff)
        self.spliney = spline(points,self.ycoeff.coeff)
        return(self.splinex.evalfull(x))
    def y(self,y):
        return(self.spliney.evalfull(y))
    def plotbasis(self):
        print("x basis:")
        self.splinex.basisplot()
        print("ybasis:")
        self.spliney.basisplot()
    def plotinter(self,pts = 1000):
        t0 = time.time()
        dt = 1/pts
        t = arange(self.basis[0], self.basis[-1] + dt, dt)
        plt.plot(self.x(t),self.y(t))
        plt.scatter(self.interx,self.intery,color="red")
        print(time.time()-t0)











points = [0,0,0,0, 0.5, 0.9, 0.94, 0.97, 1,1,1,1]

points = [0,1,2,3,4]
xy = [[2, 5], [4, 8], [4, 5], [2, 7], [3, 5], [6, 7], [4, 5], [2, 7]]
x = list(x[0] for x in xy)
y = list(y[1] for y in xy)



a = spline(points)
an = a.N


a.basisplot()
#d = d_iterative(points, 4, coeffs,a.xi)
#c = matrixequation(x,points).Avector()
#d = matrixequation(y,points).Avector()

#t0 = time.time()
#interpol= interpolation(points,x,y)

#interpol.plotbasis()
#interpol.plotinter()


#dt = 0.001
#t = arange(0, 1 + dt, dt)
#    
#ux = interpol.x(t)
#uy = interpol.y(t)
#plt.scatter(x,y,color="red")
#plt.plot(ux, uy)
#print(time.time()-t0)

#print(a)

#plt.xlim(0,6)
#plt.ylim(0,6)

#print(ux,uy)
#print(len(x),len(y))

#plt.plot(x,y)

#a = interpolate(x,y,)
#print(a.N)
#t = np.linspace(0.1,1,50)
#ty = a.evalfull(t)