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
    def __init__(self,deg,k):
        self.k = k
        self.deg = deg

    def __call__(self,x):
        return(sp.special.binom(self.deg, self.k)*
               (x)**self.k*(1-x)**(self.deg-self.k))

class bernstinepoly:
    def __init__(self,deg):
        self.deg = deg
        self.basis = self.polys()

    def polys(self):
        """
        generates a list with all the bases.
        """
        basis = []
        for k in range(self.deg+1):
            basis.append(bernsteinbasis(self.deg,k))
        return(basis)

    def __call__(self,x):

        _sum = 0
        for base in self.basis:
            _sum += base(x)

        return(_sum)

    def callbase(self,x,i):
        return(self.basis[i](x))

x = np.linspace(0,1,50)

bern = bernstinepoly(4)

y = bern(x)
plt.plot(x,y)
