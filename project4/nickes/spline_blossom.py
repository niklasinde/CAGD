#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 19:03:11 2017
 
@author: niklasinde
"""
import sys

import numpy as np
from sympy.utilities.lambdify import lambdify
import sympy as sm
import matplotlib.pyplot as plt
import time as time
import scipy as sp
plt.rcParams["figure.figsize"] = [5, 5]


class Spline:
    """ Class to calculate a spline"""

    def __init__(self, knot, k=3, coeff=None):
        if len(knot) <= k + 1:
            raise SyntaxError("We need atleast k+1 number of points in the knotvector")
        self.k = k
        self.points = knot
        if not type(coeff) == np.ndarray:
            self.coeff = [1 for i in range(len(self.points) - self.k - 1)]
        else:
            self.coeff = coeff
        self.x = sm.Symbol("x")
        self.N, self.F = self.basicfunction(self.points)

    # Loading functions
    """ functions that will load when a variable is called with this class"""

    def ref(self, p):
        """creates a list with values between the x values in points
        because sympy.Symbol("x") cant select a value and we need a x-value
        in the recurrence function"""
        return [p[i] + (p[i + 1] - p[i]) / 2 for i in range(len(p) - 1)]

    def __getrelpts(self, knot):
        n = len(knot) - self.k - 1
        p = [knot[i:i + self.k + 2] for i in range(n)]
        return (p)

    def rec(self, p, i, k, xi):
        """recurrence formula for b-spline"""
        # We should lambdify this function for higer speed.
        if k == 1:
            if p[i] == p[i + 1]:
                return 0
            if p[i] <= xi <= p[i + 1]:
                return 1
            else:
                return 0
        else:
            def div(lhs, rhs):
                if rhs == 0:
                    return 0
                else:
                    return lhs / rhs

            u = (div((self.x - p[i]), (p[i + k - 1] - p[i])) * self.rec(p, i, k - 1, xi) +
                 div((p[i + k] - self.x), (p[i + k] - p[i + 1])) * self.rec(p, i + 1, k - 1, xi))
            return u

    def update_basis(self, newknot):
        """
        Updates the spline whith a new knot vector
        Returns the index of the basisfunctions that will be changed
        """
        if len(newknot) != len(self.points):
            raise ValueError("len(newpoints) must be the same as len(oldpoints)")
        N1, F1 = self.N[:], self.F[:]

        self.points = newknot
        self.N, self.F = self.basicfunction(newknot)

        changelist = []
        for i in range(len(F1)):
            if self.F[i] != F1[i]:
#                print(True)
                changelist.append(i)
        return (changelist)

    def updatepoints(self, old, new):
        newpts = self.points[:]  # to make deep copy
        newpts.append(new)
        newpts.remove(old)
        newpts = sorted(newpts)
        return newpts

    def moveknot(self, old, new):
        newpts = self.updatepoints(old, new)
        oldrelpts = self.__getrelpts(self.points)
        newrelpts = self.__getrelpts(newpts)
        changes = list()
        for i in range(len(newrelpts)):
            if oldrelpts[i] != newrelpts[i]:
                changes.append(i)
        self.__update_basisfunction(changes, newpts)
        self.points = newpts
        return changes

    def __update_basisfunction(self, listindex: list, newpts: list) -> None:
        p = self.__getrelpts(newpts)
        for i in listindex:
            n2 = []
            f2 = []
            xi = self.ref(p[i])
            for j in range(self.k + 1):
                func = self.rec(p[i], 0, self.k + 1, xi[j])
                evalfunc = lambdify(self.x, func, modules=['numpy'])
                n2.append(evalfunc)
                f2.append(func)
            self.N[i] = n2
            self.F[i] = f2

    def basicfunction(self, knot):
        n1 = []
        f1 = []
        p = self.__getrelpts(knot)
        #        print(knot)
        for i in range(len(p)):
            n2 = []
            f2 = []
            xi = self.ref(p[i])
            for j in range(self.k + 1):
                func = self.rec(p[i], 0, self.k + 1, xi[j])
                evalfunc = lambdify(self.x, func, modules=['numpy'])
                n2.append(evalfunc)
                f2.append(func)
            n1.append(n2)
            f1.append(f2)

        return n1, f1

    def basisplot(self, i=None):
        """ Plots all the basis functions
        param i: list of the basis functions we whats to plot
        """
        if type(i) != list:
            n = range(len(self.N))
        else:
            n = i
        for i in n:
            for j in range(len(self.N[i])):
                if self.N[i][j] != 0:
                    x = np.linspace(self.points[i + j], 
                                 self.points[i + j + 1], 50)
                    y = self.coeff[i] * self.N[i][j](x)
                    plt.plot(x, y)
        plt.title("Plot of the basic functions for the splines")

    def evalfull(self, X):
        if np.shape(np.array(X)) != ():
            self.test = True
            return [self.evalfull(x) for x in X]

        hot_interval = self.hotinterval(X)
        func_val = 0
        for i in range(len(self.N)):
            if hot_interval < i or i + self.k < hot_interval:
                pass
            else:
                evalfunc = self.N[i][hot_interval - i](X)
                func_val += self.coeff[i] * evalfunc
        return func_val

    def eval_vector(self, control, xlist):
        y_list = []
        for x in xlist:
            
            hi = self.hotinterval(x)
            print(hi)
            func_val = np.array([0., 0.])
            for i in range(len(self.N)):
                if hi < i or i + self.k < hi:
                    pass
                else:
                    eval_func = self.N[i][hi - i](x)
                    func_val += control[i] * np.array(eval_func)
            y_list.append(func_val)
        y_list = np.array(y_list)
        return y_list

    def eval_basis(self, x, i):
        hot_interval = self.hotinterval_(x)
        if hot_interval < i or i + self.k < hot_interval:
            return 0.
        else:
            evalfunc = self.coeff[i] * self.N[i][hot_interval - i](x)
            return evalfunc

    def evalbasisi(self, i, x):
        return self.N[i](x)

    def hotinterval(self, x):
        p = self.points
        for i in range(len(p)-1):
            if p[i] <= x and x <= p[i+1]:
                return i
        raise ValueError(x, " is not in the interval")


class MatrixEquation:
    """ Calculates A d = x  """
    def __init__(self, xy, points):
        self.xy = xy
        self.p = points
        self.N = Spline(points)

        self.xyvec = self.xyvector()
        self.e = self.elist()
        self.A = self.Avector()

        self.coeff = self.solver()

    def Avector(self):
        """creats a square tridiagonal matrix of size n"""
        n = len(self.N.N)
        A = np.zeros((n, n))
        for col in range(n):
            for row in range(n):
                A[row, col] = self.N.eval_basis(self.e[row], col)
        return A

    def elist(self):
        p = self.p
        return list((p[i] + p[i + 1] + p[i + 2]) / 3  for i in range(1, len(self.N.N) + 1))
    def xyvector(self):
        """ When i write xy its because its the same function for x and y
        :rtype: np.array
        """
        xy = np.array(self.xy)
        xy.reshape((1, len(self.xy)))
        return xy

    def solver(self):
        return sp.linalg.solve(self.A, np.array(self.xy).transpose())


class Interpolation:
    def __init__(self, points, x, y, k=3):
        self.k = k
        self.interx = x
        self.intery = y
        self.pts = points
        self.xcoeff = MatrixEquation(x, self.pts)
        self.ycoeff = MatrixEquation(y, self.pts)
        self.splinex = Spline(self.pts, k=self.k, coeff=self.xcoeff.coeff)
        self.spliney = Spline(self.pts, k=self.k, coeff=self.ycoeff.coeff)

    def x(self, x):
        return (self.splinex.evalfull(x))

    def y(self, y):
        return (self.spliney.evalfull(y))

    def plotbasis(self):
        print("xbasis:")
        self.splinex.basisplot()
        print("ybasis:")
        self.spliney.basisplot()

    def plotinter(self, pts=1000):

        dt = 1 / pts
        t = np.arange(self.pts[0], self.pts[-1] + dt, dt)
        plt.plot(self.x(t), self.y(t))
        plt.scatter(self.interx, self.intery, color="red")


class DeBoor:
    def __init__(self, knot, control, deg):
        self.control = control
        self.knot = knot
        self.deg = deg
    def __call__(self, x):
        if np.shape(np.array(x)) != ():
            return np.array(list((self.__call__(bajs) for bajs in x)))
        return self.deBoor(self.hotinterval(x), x, self.knot, self.control, self.deg)

    def deBoor(self, k, x, t, c, p):

        """
        Evaluates S(x).

        Args
        ----
        k: index of knot interval that contains x
        x: position
        t: array of knot positions, needs to be padded as described above
        c: array of control points
        p: degree of B-spline
        """
        d = [c[j + k - p] for j in range(0, p + 1)]
        for i in range(1, p + 1):
            for j in range(p, i - 1, -1):
                
                if t[j + 1 + k - i] - t[j + k - p] == 0:
                    alpha = 1
                else:
                    alpha = (x - t[j + k - p]) / (t[j + 1 + k - i] - t[j + k - p])
                d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j]
        return d[p]

    def hotinterval(self, x):
        p = self.knot
#        print(x,p[0],"fel")
        for i in range(len(p)-1):
            if p[i] <= x and x <= p[i+1]:
                return i
        raise ValueError(x, " is not in the interval")




#knot = [0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 10]
# knot2 = [0, 0, 0, 0, 1, 2, 3, 3, 5, 6, 7, 7, 7, 7]
#knot2 = [0, 0, 0, 0, 1, 2, 3, 5, 5.5, 6, 7, 8, 9, 10, 10, 10, 10]
#s = Spline(knot)
#s.basisplot()
#s2 = Spline(knot)
#s2.basisplot()
#print("incorrect", s.moveknot(4, 5.5))
# print(s.moveknot(1,9.5))
#changes = s2.update_basis(knot2)
#print("correct  ", changes)
#s.basisplot()
#plt.ylim([0, 1])
#plt.show()

#print("spline2")
#s2 = Spline(knot2)
#s2.basisplot()
#plt.ylim([0, 1])
#plt.show()
#print("done")


knot = [1,1,1,1,6/5,7/5,8/5,9/5,2,2,2,2]
control = np.array([[0.7,-0.4],[1.0,-0.4],[2.5,-1.2],
                    [3.2,-.5],[-0.2,-.5],[.5,-1.2],[2.0,-.4],[2.3,-.4]])
#                    [0.7,−0.4],[1.0,−0.4],[2.5,−1.2],
#                    [3.2,−.5],[−0.2,−.5],[.5,−1.2],[2.0,−.4],[2.3,−.4]
save = False
D = DeBoor(knot,control,2)
A = Spline(knot, k=2)
for i in range(len(knot)-1):
    X = np.linspace(knot[i],knot[i+1],30)
#    print(X[0],"x")
    Y = D(X)
#    if i = 1
#    print(Y[:,0],Y[:,1])
    plt.plot(Y[:,0],Y[:,1])
plt.scatter(control[:,0],control[:,1])


plt.grid()
if save:
   plt.savefig("interpol2.pdf")

plt.show()


knot = [0, 0, 0, 0.3, 0.5, 0.5, 0.6, 1, 1, 1]
control = np.array([[0,0],[2,3],[7,5],[9,2],[13,1],[10,1],[7,1]])

D = DeBoor(knot,control,2)
x= np.linspace(0,0,1)






