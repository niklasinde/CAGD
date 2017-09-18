#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:41:42 2017

@author: niklasinde
"""
from numpy import *
from scipy.special import binom
from matplotlib.pyplot import *
from functools import reduce
import time

class RecurBasis:
    '''Used for task 3, to mimic the
       computation of homework 1.'''
    def __init__(self, pts):
        self.pts, self.n = pts, len(pts)-1
        self.b = self.__getCurve(0, self.n)
        self.domain = [min(self.pts[:, 0]), max(self.pts[:, 0])]

    def __getCurve(self, s, v) -> 'func':
        c1 = lambda t: (self.domain[1] - t)/(self.domain[1] - self.domain[0])
        c2 = lambda t: (t - self.domain[0])/(self.domain[1] - self.domain[0])

        def b(i, k):
            if k == 0: return lambda _: self.pts[i, :]
            else: return lambda t: c1(t)*b(i, k-1)(t) + c2(t)*b(i+1, k-1)(t)

        return b(s, v)

    def intersections(self, other:'obj', nsp=100):
        values = array(list(other.b(x)
                            for x in linspace(other.domain[0],
                                              other.domain[1], nsp)))
        xaxis = values[:, 0]
        print(xaxis)
        yaxis = values[:, 1]

        tmp = array(list(self.b(x) for x in xaxis))
        xtmp = tmp[:, 0]
        ytmp = tmp[:, 1]
        for i in range(len(tmp)):
            if abs(yaxis[i] - ytmp[i]) < 0.1:
                scatter(xaxis[i], yaxis[i])

    def subdivision(self, t=0.5) -> 'obj, obj':
        pts1 = array(list(reversed(list(self.__getCurve(0, tmp)(t)
                     for tmp in reversed(range(self.n+1))))))
        pts2 = array(list(self.__getCurve(self.n - tmp, tmp)(t)
                     for tmp in reversed(range(self.n+1))))

        A1 = RecurBasis(pts1)
        A2 = RecurBasis(pts2)
        return A1, A2

    def render(self, nsp=100, colour='r', env=True) -> 'None':
        domain = self.domain
        values = array(list(self.b(x)
                            for x in linspace(domain[0], domain[1], nsp)))

        plot(values[:, 0], values[:, 1], c=colour)

        if env:
            print('env active')
            LineSegment = lambda i, j: [self.pts[i, j], self.pts[i+1, j]]
            scatter(self.pts[:, 0], self.pts[:, 1], c='k', alpha=0.5)
            plot(list(LineSegment(i, 0) for i in range(self.n)),
                 list(LineSegment(i, 1) for i in range(self.n)),
                 c='k', alpha=0.2)


