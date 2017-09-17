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
#            return x
            return (x - do[0])/(do[1]-do[0])
        def one_m_t(x):
#            return 1-x
            return (do[1]-x)/(do[1]-do[0])
        var = lambda x: (t(x)**self.i * one_m_t(x)**(self.n - self.i))
        return lambda arg: self.coeff*var(arg)

class BernPoly:
    def __init__(self, pts:'list/vector'):
        self.pts, self.n = array(pts), len(pts) - 1
        self.domain = [min(pts[:,0]),max(pts[:,0])]
        self.__lineartrans(pts)
#        print(self.transpoints)
        self.basis = self.__getBasis()
        self.B     = self.__getPolynomial()

    def __call__(self, x) -> 'matrix':
        return self.B(x)

    def __getPolynomial(self) -> 'func':
        """ Lamdyfies __getBasis """
#        print(self.transpoints)
        return lambda x: sum(self.pts[i]*base(x)
                             for i, base in enumerate(self.basis))

    def __getBasis(self) -> 'matrix':
        """ Calls the BernBasis class """
        return array(list(BernBasis(self.domain, self.n, idx) for idx in range(self.n + 1)))


    def __lineartrans(self,pts):
        do = self.domain
        def t(x):
            return (x - do[0])/(do[1]-do[0])
        for i in range(len(pts)):
            pts[i,0] = t(pts[i,0])
        self.transpoints = pts

    def render(self, nsp=100, colour='r', env=True) -> 'None':
        """ Plotting function"""
        domain = self.domain
#        print(domain)
        sample = linspace(domain[0], domain[1], nsp)
        values = array(list(self.B(x)
                       for x in linspace(domain[0], domain[1], nsp)))

        plot(values[:, 0], values[:, 1], c=colour, alpha=0.5)

        if env:
            LineSegment = lambda i, j: [self.pts[i, j], self.pts[i+1, j]]
            scatter(self.pts[:, 0], self.pts[:, 1], c='k', alpha=0.5)
            plot(list(LineSegment(i, 0) for i in range(self.n)),
                 list(LineSegment(i, 1) for i in range(self.n)),
                 c='k', alpha=0.2)

    def renderBasisFuncs(self, pt = False) -> 'None':
        tmp = self.basis.copy()
        sample = linspace(0, 1, 100)
        for i, p in enumerate(tmp):
            plot(sample, vectorize(p)(sample))

        xlim(0,1)
        ylim(-0.2, 1.2)

        if pt: print(sum(base(pt) for base in tmp))

    def subdivision(self, t=0.5):
        def b(i, k):
            if k == 0: return lambda _: self.pts[i, :]
            else: return lambda t: (1-t)*b(i, k-1)(t) + t*b(i+1, k-1)(t)

        pts1 = array(list(reversed(list(b(0, tmp)(t)
                     for tmp in reversed(range(self.n+1))))))
        pts2 = array(list(b(self.n - tmp, tmp)(t)
                     for tmp in reversed(range(self.n+1))))

        scatter(pts1[:, 0], pts1[:, 1], c='r')
        scatter(pts2[:, 0], pts2[:, 1], c='b')

        A1 = RecurBasis(pts1)
        A1.render(colour = 'r')
        A2 = RecurBasis(pts2)
        A2.render(colour = 'b')
