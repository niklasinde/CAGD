#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 12:18:24 2017

@author: niklasinde
"""

import scipy as sp
import numpy as np
from sympy.utilities.lambdify import lambdify
import sympy as sm
import matplotlib.pyplot as plt


class bernsteinbasis:
    def __init__(self, deg, k):
        self.k = k
        self.deg = deg

    def __call__(self, x):
        return (sp.special.binom(self.deg, self.k)*
                (x)**self.k*(1 - x)**(self.deg - self.k))

    def __call__(self,x):
        return(sp.special.binom(self.deg, self.k) * (x)**self.k * (1-x)**(self.deg-self.k))
        
        
    def __repr__(self):
        x = sm.Symbol("X")
        return(str(sp.special.binom(self.deg, self.k) * (x)**self.k * (1-x)**(self.deg-self.k)))
        

class bernsteinpoly:
    def __init__(self, deg):
        self.deg = deg
        self.basis = self.polys()

    def __call__(self,x):
        return sum(self.pts[i]*base(x)
                   for i, base in enumerate(self.basis))

        self.coeff = np.ones((1,self.deg+1))

    def polys(self):
        """
        generates a list with all the bases.
        """

        return list(bernsteinbasis(self.deg, k) for k in range(self.deg + 1))

    def callbase(self,x,i):
        return(self.basis[i](x))

        basis = [bernsteinbasis(self.deg,k) for k in range(self.deg+1)]
        return basis
            
    def __call__(self,x):
        
        return sum(base(x) for base in basis)
        for base in self.basis:
            _sum += base(x)
            
        return(_sum)
        
    def callbase(self,x,i):
        return(self.basis[i](x))


    


x = np.linspace(0,1,50)

bern = bernstinepoly(4)


y = bern(x)
        
        
for i in range(4):
    plt.plot(x,bern.callbase(x,i))
>>>>>>> Stashed changes
plt.plot(x,y)
