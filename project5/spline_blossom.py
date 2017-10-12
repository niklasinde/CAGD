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
plt.rcParams["figure.figsize"] = [9, 5]

          
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
#        print(k,t[k])
#        print(c[p+1 + k - p])
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
        for i in range(len(p)-1):
            if p[i] <= x and x <= p[i+1]:
                return i
        raise ValueError(x, " is not in the interval")
    def render(self):
        for i in range(len(self.knot)-1):
            X = np.linspace(self.knot[i],self.knot[i+1],30)

            Y = self.__call__(X)
            plt.plot(Y[:,0],Y[:,1])
        plt.scatter(control[:,0],control[:,1])
        plt.show()

        
class Spline:
    """ Class to calculate a spline"""

    def __init__(self, knot, k=3, coeff=None): 
        """
        Knot is the knot vector
        Coefficients are inputed if we have 1d coeffcents else use the eval_vector function.
        """
        if len(knot) <= k + 1:
            raise SyntaxError("We need atleast k+1 number of points in the knotvector , now you have "+str(len(knot
                                                                                                        ))+" and you need "+str(k+1))
        if sorted(knot) != knot:
            raise Exception("knot needs to be sorted.")
        self.k = k # degree of basis spline
        self.points = knot
        if not type(coeff) == np.ndarray:
            self.coeff = [1 for i in range(len(self.points) - self.k - 1)]
        else:
            self.coeff = coeff
        
        self.x = sm.Symbol("x")
        self.N, self.F = self.basisfunction(self.points)

    # Loading functions
    """ functions that will load when a variable is called with this class"""

    def ref(self, p):
        """creates a list with values between the x values in points
        because sympy.Symbol("x") cant select a value and we need a x-value
        in the recurrence function"""
        return [p[i] + (p[i + 1] - p[i]) / 2 for i in range(len(p) - 1)]

    def __getrelpts(self, knot):
        """
        Relative points are used for calculating the b-spline
        The B-spline, B_i(x) have the knots p_i
        By using this we can have i = 0 in the reccurent formula.
        """
        n = len(knot) - self.k - 1
        p = [knot[i:i + self.k + 2] for i in range(n)] # len(p[0]) = self.k + 2
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
    def basisfunction(self, knot):
        """
        Creates the basisfunctions
        
        Each basis function consists of its own functions.
        They are stored in n2 and then appended to n1.
        
        This functions then returns n1 and f1. n1 is the lambdified version of f1
        but its not used for computing only if we wich to see the basis equations.
        
        """
        n1 = list()
        f1 = list()
        p = self.__getrelpts(knot)
        for i in range(len(p)):
            n2 = list()
            f2 = list()
            xi = self.ref(p[i])
            for j in range(self.k + 1):
                func = self.rec(p[i], 0, self.k + 1, xi[j])
                evalfunc = lambdify(self.x, func, modules=['numpy'])
                n2.append(evalfunc)
                f2.append(func)
            n1.append(n2)
            f1.append(f2)
        return n1, f1

    def __updatepoints(self, old, new):
        """
        Updates the knot vector
        """
        newpts = self.points[:]  # to make deep copy
        newpts.append(new)
        newpts.remove(old)
        newpts = sorted(newpts)
        return newpts

    def moveknot(self, old, new):
        """
        Because of local support we do not have to change all the basis functions.
        Only the basis function whos not have changed.
        Instead of commparing the relative points we could find the change
        by knowing the index of the new and old points.
        """
        newpts = self.__updatepoints(old, new)
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
        """
        Only changes the basisfunction we want
        """
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

    def evalfull(self, X):
        """
        Evaluates the sum{c_iB_i(x)}
        """
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
        """
        eval_full but the coefficents 2d. (controlpoints)
        """
        if len(self.N) != len(control):
            raise Exception("You need "+ str(len(self.N)) + " controlpoints. Now you have "+ str(len(control)))
        y_list = []

        for x in xlist:
            hi = self.hotinterval(x)
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

    def render_vector(self, control):
        """
        plots eval_vector with different colors in each segment
        """
        if len(self.N) != len(control):
            raise Exception("You need "+ str(len(self.N)) + " controlpoints. Now you have " + str(len(control)))
        FullSupport = lambda knot, t, deg: knot[deg] <= t <= knot[-deg-1]
        for i in range(len(knot)-1):

            if not (FullSupport(knot,knot[i], 3) and FullSupport(knot, knot[i+1], 3)):
                X = np.linspace(knot[i],knot[i+1],30)
                Y = self.eval_vector(control,X)
                plt.plot(Y[:, 0],Y[:, 1], color = "black",linewidth=2)
            else:
                X = np.linspace(knot[i],knot[i+1],30)
                Y = self.eval_vector(control,X)
                plt.plot(Y[:, 0],Y[:, 1],linewidth=2)
        plt.scatter(control[:, 0],control[:, 1])
#        plt.show()
        
    
    def eval_basis(self, x, i):
        hot_interval = self.hotinterval(x)
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
    def check_full_support(self):
        FullSupport = lambda knot, t, deg: knot[deg] <= t <= knot[-deg-1]
        
        for i in range(len(self.points)):
            if FullSupport(self.points,self.points[i],self.k):
                pass
            else:
                print("You dont have full support")
                return(False)
        return True
      




# Task 1
#knot = [0,0,0,0,0,1,2,3,4,5,6,7,8,8,8,8,8]
#control = np.array([[4,0], [2,1], [2,1], [8,10], [4,2],
#                    [10,0],[5,-1],[4,-7],[2,-4],[2,1],[-1,2],[10,10]])
#control2 = np.array([[4,0], [2,1], [2,1], [8,10], [4,2],
#                    [5,0],[5,-1],[4,-7],[2,-4],[2,1],[-1,2],[10,10]])
#A = Spline(knot, k=4, coeff=control[:,0])
#B = Spline(knot, k=4, coeff=control2[:,0])
#
#O = Spline(knot, k=4)
##O.basisplot()
##plt.savefig("bajs.pdf")
#
#A.basisplot()
#plt.savefig("Basis1.pdf")
#plt.show()
##
#B.basisplot()
#plt.savefig("Basis2.pdf")
#plt.show()
#A.render_vector(control)
#plt.savefig("Interpol1.pdf")
#plt.show()
#A.render_vector(control2)
#plt.savefig("Interpol2.pdf")
#plt.show()
        
    
# Task 5
#knot = [0, 1/11, 2/11, 3/11, 4/11, 5/11, 6/11, 7/11, 8/11, 9/11, 10/11, 1]
#control = np.array([ [0,0], [3,2], [9,-2], [7,-5], [1,-3], [1,-1], [3,1], [9,-1]]) 
##control2 = np.array([ [0,0], [3,2], [9,-2], [7,-5], [1,-3], [1,-1], [3,1], [0,0]]) 
##control3 = np.array([ [0,0], [3,2], [9,-2], [7,-5], [1,-3], [1,-1], [0,0],[9,1]]) 
##control4 = np.array([ [0,0], [3,2], [9,-2], [7,-5], [1,-3], [0,0], [3,1],[9,1]]) 
##
#A = Spline(knot,k=3,coeff = control[:,0])
#A.render_vector(control)
##plt.savefig("Task51.pdf")
#plt.show()
#A.basisplot()
##plt.savefig("Task51basis.pdf")
#plt.show()
#
#B = Spline(knot,k=3)
#B.render_vector(control2)
#plt.savefig("Task52.pdf")
#plt.show()
#
#B = Spline(knot,k=3)
#B.render_vector(control3)
#plt.savefig("Task53.pdf")
#plt.show()
#
#B = Spline(knot,k=3,coeff = control4[:,0])
#B.render_vector(control4)
#plt.savefig("Task54.pdf")
#plt.show()


class MatrixEquation:
    """ Calculates A d = x  """
    def __init__(self, xy, points):
        self.xy = xy
        self.p = points
        self.N = Spline(points)
        print(len(self.N.N),"lenN.N")
        print(len(self.p),"lenpts")
        self.xyvec = self.xyvector()
        self.e = self.elist()
        print(len(self.e),"len of e")
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
        return list((p[i] + p[i + 1] + p[i + 2]) / 3  for i in range(len(self.N.N)))
    def xyvector(self):
        """ When i write xy its because its the same function for x and y
        :rtype: np.array
        """
        xy = np.array([self.xy])
        xy.reshape((1, len(self.xy)))
        return xy

    def solver(self):
        print(self.A.shape,self.xyvec.shape)
        return np.linalg.solve(self.A, self.xyvec.T)


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

    def plotinter(self, pts=100):
        dt = 1 / pts
        t = np.arange(self.pts[0], self.pts[-1] + dt, dt)
        plt.plot(self.x(t), self.y(t))
        plt.scatter(self.interx, self.intery, color="red")




xy = [[0,0],[6,10],[7,10.2],[9,8]]
x = [x[0] for x in xy]
y = [x[1] for x in xy]

def unipoints(a,b,s, k):
    l = [a] * k
    for k in range(1,s-1):
        tk = a+k*(b-a)/s
        l.append(tk)
    for i in range(k):
        l.append(b)
    return(l)
k = 3
s = 5
uni = unipoints(0,1,s,k+1)
print(x,"x")
print(y,"y")
print(unipoints(0,1,s,k+1))

A = Interpolation([0,0,0,1,2,3,6,6],x,y,k=k)
A.plotinter()


