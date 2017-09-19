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
        self.pts, self.n = array(pts, copy =True), len(pts) - 1
        self.domain = [min(pts[:,0]),max(pts[:,0])]
        self.__lineartrans()
        self.basis = self.__getBasis()
        self.B     = self.__getPolynomial()
        self.testlist()

       
    def testlist(self):
        self.l = []

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


    def __lineartrans(self):
        """ linear transformation of points to be able to have any domain"""
        do = self.domain
        self.transpoints = copy(self.pts)
        def t(x):
            return (x - do[0])/(do[1]-do[0])
        for i in range(len(self.transpoints)):
            self.transpoints[i,0] = t(self.transpoints[i,0])

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
        self.render(env=env, nsp=nsp, color=color)
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
    def subdiv(self, ti=0.5):
        do = self.domain
        if not do[0] <= ti <= do[1]:
            raise ValueError(ti,do)
        def f(x):
            return (x - do[0])/(do[1]-do[0])
        t = f(ti)
        def inv(x):
            return(x*(do[1]-do[0])+do[0])
        def b(i, k):
            if k == 0: return lambda _: self.pts[i, :]
            else: return lambda t: (1-t)*b(i, k-1)(t) + t*b(i+1, k-1)(t)

        self.pts1 = array(list(reversed(list(b(0, i)(t)
                     for i in reversed(range(self.n+1))))))
        self.pts2 = array(list(b(self.n - i, i)(t)
                     for i in reversed(range(self.n+1))))

        self.A1 = BernPoly(self.pts1)
        self.A2 = BernPoly(self.pts2)
        return(self.A1,self.A2)
    def bezint(self,other):
        threshold=19**-10
        """
        Other == Line
        Does bbox(B1) intersect bbox(B2)?
        No: Return false. 
        Yes: Continue.
        Is area(bbox(B1)) + area(bbox(B2)) < threshold? done
        Yes: Return true.
        No: Continue.
        Split B1 into B1a and B1b at t = 0.5
        Split B2 into B2a and B2b at t = 0.5
        Return bezInt(B1a, B2a) || bezInt(B1a, B2b) || bezInt(B1b, B2a) || bezInt(B1b, B2b).
        """
#        other.pts2 = np.array(list(other(x) for x in linspace(self.domain[0],self.domain[1],100)))
        B1, B2 = ch.convexhull(self.pts), ch.convexhull(other.pts)
        if abs(self.domain[0]-self.domain[1])<threshold and B1.intersect(B2):
            other.l.append(self)
        
        elif B1.intersect(B2): 
#            self.render()
#            other.plot()
            B1a,B1b = self.subdiv(self.domain[0]+(abs(self.domain[1]-self.domain[0]))/2)
            return(B1a.bezint(other), B1b.bezint(other))
        else:
            return False

    def bezintwrapp(self,other,plot=True):
        A = self.bezint(other)
        if plot:
            self.render()
            other.render()
            for i in B.l:
                scatter(i(i.domain)[0],i(i.domain)[1],color ="g")
#                print(i.domain)
        x = array(list(i(i.domain[0]) for i in other.l))
        return(x)
    

        
        
        
        
        
    
    
points1 = array([[ -1   ,  0  ],[ 0 ,  1   ],[ 1, -2   ], [ 2 , 0  ]])
points2 = (array([[0.,0.],[9.,-4.],[7.,5.],[2.,-4.]]))
points3 = array([[4,5],[6,-4]])
A = BernPoly(points2)

B = BernPoly(points3)
#print(type(A))
#A.render()
#B.render()
#a,b = A.subdiv(0.5)
#a.render()
#b.render()
A.bezintwrapp(B)

#bajs = B.l
#x = [x(x.domain[0])[0] for x in bajs]
#y = [x(x.domain[0])[1] for x in bajs]
#A.render()
#B.render()
#scatter(x,y,color="red")
show()
#print(x)
#print(y)
#print(A.bezint(B),"result")
#Ahull=ch.convexhull(A.pts)
#Bhull=ch.convexhull(B.pts)

#print(Ahull.intersect(Bhull),"test")

#Bhull.render()
#Ahull.render()
#show()









