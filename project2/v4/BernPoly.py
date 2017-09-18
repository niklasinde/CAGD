#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:38:45 2017

@author: niklasinde
"""

from numpy import *
from scipy.special import binom
from matplotlib.pyplot import *
from functools import reduce
import time
import convexhull as ch

def timeit(method):
    def timed(*args, **kw):
        ts     = time.time()
        result = method(*args, **kw)
        te     = time.time()

        print('%r (%r, %r) %2.2f sec' % (method.__name__, args, kw, te-ts))
        return result
    return timed


class BernBasis:
    def __init__(self, domain ,degree:'int', iteration:'int'):
        
        self.n, self.i, self.domain = degree, iteration, domain
        assert self.n >= self.i, 'passed max iter.'
        self.coeff = binom(self.n, self.i)
        self.evaluate = self.__getB()


    def __call__(self, arg:'float') -> 'float':
        return self.evaluate(arg)
    

    def __getB(self) -> 'float':
        do = self.domain
        def t(x):
            return (x - do[0])/(do[1]-do[0])
        def one_m_t(x): #one minus t
            return (do[1]-x)/(do[1]-do[0])
        var = lambda x: (t(x)**self.i * one_m_t(x)**(self.n - self.i))
        return lambda arg: self.coeff*var(arg)

class BernPoly:
    def __init__(self, pts:'list/vector'):
#        print(len(pts),"len(points)")
            
        self.pts, self.n = array(pts), len(pts) - 1
        self.domain = [min(pts[:,0]),max(pts[:,0])]
        self.__lineartrans(pts)
        self.basis = self.__getBasis()
        self.B     = self.__getPolynomial()

    def __call__(self, x) -> 'matrix':
        """
        Call function: Returns the x and y coordinates.
        """
        return self.B(x)

    def __getPolynomial(self) -> 'func':
        """ Lamdyfies __getBasis and adds then coefficients(points)"""
        return lambda x: sum(self.pts[i]*base(x)
                             for i, base in enumerate(self.basis))

    def __getBasis(self) -> 'matrix':
        """ Calls the BernBasis class """
        return array(list(BernBasis(self.domain, self.n, idx) for idx in range(self.n + 1)))


    def __lineartrans(self,pts):
        """ linear transformation of points to be able to have any domain"""
        do = self.domain
        def t(x):
            return (x - do[0])/(do[1]-do[0])
        for i in range(len(pts)):
            pts[i,0] = t(pts[i,0])
        self.transpoints = pts

    def render(self, env=True, nsp=100, color='r') -> 'None':
        """ Plotting function"""
        domain = self.domain
        sample = linspace(domain[0], domain[1], nsp)
        values = array(list(self.B(x)
                       for x in linspace(domain[0], domain[1], nsp)))
        plot(values[:, 0], values[:, 1], c=color, alpha=0.5)
        if env:
            LineSegment = lambda i, j: [self.pts[i, j], self.pts[i+1, j]]
            scatter(self.pts[:, 0], self.pts[:, 1], c='k', alpha=0.5)
            plot(list(LineSegment(i, 0) for i in range(self.n)),
                 list(LineSegment(i, 1) for i in range(self.n)),
                 c='k', alpha=0.2)
    def plot(self,env=True,nsp=100,color='r'):
        self.render(env=env, nsp=nsp, colour=color)
        show()
    def renderConvexhull(self):
        con = ch.convexhull(self.pts)
        con.render()
    def plotConvexhull(self):
        con = ch.convexhull(self.pts)
        con.plot()
    def renderBasisFuncs(self) -> 'None':
        tmp = self.basis.copy()
        sample = linspace(0, 1, 100)
        for i, p in enumerate(tmp):
            plot(sample, vectorize(p)(sample))
        xlim(self.domain[0],self.domain[1])
        ylim(-0.2, 1.2)
    def subdiv(self, t, overwrite =None):
        do = self.domain
        def f(x):
            return (x - do[0])/(do[1]-do[0])
        t = f(t)
        def b(i, k):
            if k == 0: return lambda _: self.pts[i, :]
            else: return lambda t: (1-t)*b(i, k-1)(t) + t*b(i+1, k-1)(t)

        pts1 = array(list(reversed(list(b(0, i)(t)
                     for i in reversed(range(self.n+1))))))
        pts2 = array(list(b(self.n - i, i)(t)
                     for i in reversed(range(self.n+1))))

        A1 = BernPoly(pts1)
        A2 = BernPoly(pts2)
        if overwrite != None:
            Return
        else:
            A1 = BernPoly(pts1)
            A2 = BernPoly(pts2)
            return(A1,A2)
    
    def intersect(self,other,threshold=1):
        """
        Does bbox(B1) intersect bbox(B2)?
        No: Return false. 
        Yes: Continue.
        Is area(bbox(B1)) + area(bbox(B2)) < threshold?
        Yes: Return true.
        No: Continue.
        Split B1 into B1a and B1b at t = 0.5
        Split B2 into B2a and B2b at t = 0.5
        Return bezInt(B1a, B2a) || bezInt(B1a, B2b) || bezInt(B1b, B2a) || bezInt(B1b, B2b).
        """
        a, b = ch.convexhull(self.pts), ch.convexhull(other.pts)
        for i in range(100):    
            if a.area() + b.area() < threshold: return(a.domain,b.domain)
            elif a.intersect(b): 
            

#
#                
#        for i in range(len(otherx)):
#            if domain[0] <= other(otherx[i])[0]<= domain[1]:
#                xo.append(otherx[i])
#        print(xs,"\n",xo)
   
        
        
        
        
        
    
    
pts1 = array([[ -1   ,  0  ],[ 0 ,  1   ],[ 1, -2   ], [ 2 , 0  ]])
pts2 = array([[0,0],[9,-4],[7,5],[2,-4]])
pts3 = array([[4,5],[6,-4]])
A = BernPoly(pts2)
B = BernPoly(pts3)


a.b = A.subdiv(0.5)



#a = A.subdiv(1.7177033492822966)[1]
#b = a.subdiv(5.5071770334928232)[0]

a.render(color="black")
print(a.domain)
#b.render(color = "green")
B.render()
















