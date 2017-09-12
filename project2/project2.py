#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 12:18:24 2017

@author: niklasinde
"""

import scipy as sp

from sympy.utilities.lambdify import lambdify
import sympy as sm


class bernsteinbasis:
    def __init__(self,deg,k):
        
        if deg > 0:
            self.k = k
            self.deg = deg
            self.multiplebasis = False
        
        self.x = sm.Symbol("x")
        self.basefuncsm = self.polybase()
        self.basefunceval = lambdify(self.x, self.basefuncsm, modules=['numpy'])

    def polybase(self):
        Bvn = sp.special.binom(self.deg, k) * (self.x)**k * (1-self.x)**(self.deg-self.k)
        
    def __call__(self,t):
        return(self.basefunceval(t))
        
    def __add__(self,other):
        
        new_class = bernsteinbasis(other.deg,-1)
        
        new_ca
        
        add = self.basefuncsm + other.basefuncsm
        return(add)
        




class bernstinepoly:
    def __init__(self,deg):
        
        self.deg = deg
        self.polylist = self.polys()
        
    def polys(self):
        """
        generates a list with all the bases.
        """
        polylist = []
        for i in range(len(deg)+1):
            polylist.append(bernsteinbasis(self.deg),i)
            
            
    def __call__(self,t):
        pass        