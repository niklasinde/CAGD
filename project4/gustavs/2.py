#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 16:46:28 2017

@author: niklasinde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 09:19:19 2017

@author: filip
"""

from  scipy import *
from  pylab import *
from functools import reduce
from sympy.utilities.lambdify import lambdify
import sympy as sm

class B_spline:
    def __init__(self, p, t, coeffs):
        
        len_splines = len(t) - 1 + p
        self.p = p
        self.coeffs = coeffs
        self.x = sm.Symbol("x")
        
        # Padd t-array
        self.t = t[:] # Copy the list in order to not destroy it by padding
        self.t[0 : 0] = [self.t[0]] * p
        self.t[-1 : -1] = [self.t[-1]] * p
        
        self.splines = []
        for i in range(len_splines):
            self.splines.append(B(self.p, self.t, i))
            
    def __call__(self, x):
        y = []
        for x_i in x:
            spline_vals = list(map(lambda B, c: c * B(x_i), self.splines, self.coeffs))
            sum_ = reduce(lambda v1, v2: v1 + v2, spline_vals)
            y.append(sum_)
        
        return y

class B:
    def __init__(self, p, t, i):
        self.p = p
        self.t = t
        self.i = i
    
    def __call__(self, x):
        return self.evaluate(self.i, self.p, self.t, x)
    
    def evaluate(self, i, p, t, x):
        if p == 0:
            return 1 * (t[i] <= x) * (x < t[i + 1])
        else:
            def div(lhs, rhs):
                return 0 if rhs == 0 else lhs / rhs
                    
            term1 = div(x - t[i], t[i + p] - t[i]) * self.evaluate(i, p - 1, t, x)
            term2 = div(t[i + p + 1] - x,  t[i + p + 1] - t[i + 1]) * self.evaluate(i + 1, p - 1, t, x)
            return term1 + term2

class ControlPoint():
    def __init__(self, knot_min, knot_max, coeff):
        self.knot_min = knot_min
        self.knot_max = knot_max
        self.coeff = coeff
        
    def __repr__(self):
        return '\nmin: {}, max {}, coeff: {}'.format(self.start, self.end, self.coeff)
    
class D():
    # k = knot
    def __init__(self, val, knots):
        self.val = val
        self.knots = knots
        
    def __repr__(self):
        if not (len(self.knots) == 0):
            return '\nval: {}, k_min: {}, k_max {}'.format(self.val, self.knots[0], self.knots[-1])
        else:
            return '\nfinal_value: {}'.format(self.val)
    
    @classmethod
    def fromOtherD(cls, d1, d2, u):
        
        a = D.alpha(d1.knots, d2.knots, u)
        val = a * d1.val + (1 - a) * d2.val
        
        knots = d1.knots[1:] # Remove the first knot from d1 and take that as knot vector
        
        return cls(val, knots)
    
    @classmethod
    def alpha(cls, d1_knots, d2_knots, u):
        k_min = d1_knots[0]
        k_max = d2_knots[-1]
        
        return (k_max - u) / (k_max - k_min)

# cps = control points
# u = eval point
# n = order of spline - NOTE: only algorithm for n = 4 is provided
def d_iterative(t, n, coeffs, u):
    
    p = n - 1
    t_ = t[:] # Copy to not modify input args
    # Padd t-array
    t_[0 : 0] = [t_[0]] * p
    t_[-1 : -1] = [t_[-1]] * p
    
    d = list(map(lambda x: d_element_wise_iteration(coeffs, t_, x, n).val, u))
    
    return d
    
def d_element_wise_iteration(c, t, u, n):
    
    if not (t[0] <= u <= t[-1]):
        return D(0, [])
    
    d_vals = extract_cps(c, t, u)
    
    rows = cols = n
    d = zeros((rows, cols), dtype = object) # It is possible to save memory here
    d[:, 0] = d_vals # Replace first row with known coeffs
    
    # Iteratively create a matrix with the values for the next column in d, with the first column initialized
    for i in range(cols - 1):
        for j in range(rows - 1 - i): # For each row which contains a value, except the last (which is only used once)
            d[j, i + 1] = D.fromOtherD(d[j, i], d[j + 1, i], u);
            
    return d[0, cols - 1]

# From all control points, extract the coeffs needed for calculation of d = s(u)
def extract_cps(c, t, u):
    
    index_u_i = get_index_u_i(t, u)
    index_u_iminus2 = index_u_i - 3 # Note: yes, 3, not 2, since that is the starting point of the spline
    
    d_vals = []
    for i in range(index_u_iminus2, index_u_iminus2 + 4, 1):
        d_vals.append(D(c[i], t[i + 1 : i + 4]))
    
    return d_vals

# Binary search to find index
# u = value between u_i and  u_i1+1
def get_index_u_i(t, u):
    
    t = t[3:-3]
    first = 0
    last = len(t) - 1
    mid = last // 2
    if t[mid] == t[-1]:
            return mid + 3
    
    u_left = t[mid]
    u_right = t[mid + 1]
    
    while (not (u_left <= u <= u_right)): # Note: <= on both sides, since nodes might be coincident
    
        if (u > u_left):
            first = mid + 1
        else:
            last = mid
            
        mid = (last + first) // 2
        
        if t[mid] == t[-1]: # The last point corresponding to Ui+1.start is missing (t[-1]), and we must ensure that u_right doesnt go out of bounds
            return mid + 3
        
        u_left = t[mid]
        u_right = t[mid + 1]
    
    return mid + 3

def view_b_spline_from_poly(p, t, coeffs, x):
    
    spline = B_spline(p, t, coeffs)
    y = spline(x)
    
    plot(x, y)
    
def view_b_spline_from_blossom(p, t, coeffs, x):
    
    n = p + 1
    d = d_iterative(t, n, coeffs, x)
    
    plot(x, d)
    
def exper3():
    
    p = 3
    t_start = 0
    t_end = 2
    
    t = list(range(t_start, t_end + 1, 1))
    x = arange(t_start - 1, t_end + 1, 0.01)
    
    len_splines = len(t) - 1 + p
    coeffs = [5] * len_splines
    
    coeffs[0] = 5
    coeffs[2] = 9
    coeffs[3] = 13
    print(coeffs,t)
    
    view_b_spline_from_poly(p, t, coeffs, x)
    view_b_spline_from_blossom(p, t, coeffs, x)
    
exper3()






