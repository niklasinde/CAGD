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
from mpl_toolkits.mplot3d import axes3d

plt.rcParams["figure.figsize"] = [9, 5]
np.set_printoptions(linewidth = 1000)


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
            raise Exception("We need atleast k+1 number of points in the knotvector , now you have "+str(len(knot
                                                                                                        ))+" and you need "+str(k+1))
        if sorted(knot) != knot:
            raise Exception("knot needs to be sorted.")

        self.k = k # degree of basis spline
        self.points = knot

        self.x = sm.Symbol("x")


        self.N, self.F = self.basisfunction(self.points)
        if not type(coeff) == np.ndarray:
            self.coeff = [1 for i in range(len(self.N))]
        else:
            self.coeff = coeff

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
                if self.F[i][j] != 0:
#                    print(self.F[i][j],"hej")
                    x = np.linspace(self.points[i + j],
                                    self.points[i + j + 1], 50)
                    y = self.N[i][j](x)
                    plt.plot(x,y)

                else:
                    x = np.linspace(self.points[i + j],
                                    self.points[i + j + 1], 50)
#                    print(coeff[i], self.N[i][j](x),"y")

    def evalfull(self, X):
        """
        Evaluates the sum{c_iB_i(x)}
        """
        if np.shape(np.array(X)) != ():
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

    def eval_vector(self, control, x):

        """
        eval_full but the coefficents 2d. (controlpoints)
        """
        if len(self.N) != len(control):
            raise Exception("You need " + str(len(self.N)) +
                            " controlpoints. Now you have "
                            + str(len(control)))

        if np.shape(np.array(x)) != ():
            return np.array([self.eval_vector(control, cx_i) for x_i in x])

        hi = self.hotinterval(x)
        func_val = np.zeros(control.shape[1])
        for i in range(len(self.N)):
            if hi < i or i + self.k < hi:
                pass
            else:
                eval_func = self.N[i][hi - i](x)
                func_val += control[i] * np.array(eval_func)
        return func_val

    def eval_nurbs(self, control, weights, x):

        if len(self.N) != len(control) or len(self.N) != len(weights):
            raise Exception("You have {} basisfunction. Now you have {} controlpoints and {} wights".format(len(self.N),len(control), len(weights)) )
        y_list = []
        if np.shape(np.array(x)) != ():
            return np.array([self.eval_nurbs(control,weights, x_i) for x_i in x])

        hi = self.hotinterval(x)
        teller = np.zeros(control.shape[1])
        denominator = 0
        for i in range(len(self.N)):
            if hi < i or i + self.k < hi:
                pass
            else:
                eval_func = self.N[i][hi - i](x)

                teller += control[i] * np.array(eval_func) * np.array(weights[i])
                denominator += np.array(eval_func) * np.array(weights[i])
        if denominator.any == 0:
            return (teller)
        return teller / denominator
    # Eval functions

    def eval_basis(self, x, i):
        hot_interval = self.hotinterval(x)
        if hot_interval < i or i + self.k < hot_interval:
            return 0.
        else:
            evalfunc = self.coeff[i] * self.N[i][hot_interval - i](x)
            return evalfunc

    def hotinterval(self, x):

        ### Make a bineary search of hot interval.
        p = self.points
        for i in range(len(p)-1):
            if p[i] <= x and x <= p[i+1]:
                return i
        raise ValueError(str(x) + " is not in the interval")

    def check_full_support(self):
        def fullSupport(knot, t, deg): return knot[deg] <= t <= knot[-deg-1]

        for i in range(len(self.points)):
            if fullSupport(self.points, self.points[i], self.k):
                pass
            else:
                print("You dont have full support")
                return(False)
        return True

    def plot2D(self, control=None, interpoints=None, pts=30):

        """
        plots eval_vector with different colors in each segment
        """
        if np.array(control).shape == ():
            control = self.coeff
        if len(self.N) != len(control):
            raise Exception("You need " + str(len(self.N)) + "controlpoints. Now you have " + str(len(control)))
        knot = self.points

        def fullSupport(knot, t, deg): return knot[deg] <= t <= knot[-deg-1]

        for i in range(len(self.points)-1):

            if not (fullSupport(knot, knot[i], self.k) and fullSupport(knot, knot[i+1], self.k)):
                X = np.linspace(knot[i], knot[i+1], pts)
                Y = self.eval_vector(control, X)
                plt.plot(Y[:, 0], Y[:, 1], color="black", linewidth=2)
            else:
                X = np.linspace(knot[i], knot[i+1], pts)
                Y = self.eval_vector(control, X)
                plt.plot(Y[:, 0], Y[:, 1], linewidth=2)
        if np.array(interpoints).shape != ():
            plt.scatter(control[:, 0],control[:, 1], color="red")
            plt.scatter(interpoints[:,0],interpoints[:, 1], color="black")
        else:
            plt.scatter(control[:, 0],control[:, 1])
        plt.show()

    def plotnurbs(self, control, weights, pts=30, show=False):

        """
        plots eval_vector with different colors in each segment
        """
        if np.array(control).shape == ():
            control = self.coeff
        if len(self.N) != len(control) != len(weights):
            raise Exception("You need {} controlpoints. Now you have {} controlpoints and {} wights".format(len(self.N),len(control), len(weights)) )
        knot = self.points

        def fullSupport(knot, t, deg): return knot[deg] <= t <= knot[-deg-1]

        for i in range(len(self.points)-1):
            X = np.linspace(knot[i], knot[i+1], pts)
            Y = self.eval_nurbs(control, weights, X)
            if not (fullSupport(knot, knot[i], self.k) and
                    fullSupport(knot, knot[i+1], self.k)):
                plt.plot(Y[:, 0], Y[:, 1], color="black", linewidth=2)
            else:
                plt.plot(Y[:, 0], Y[:, 1], linewidth=2)
        else:
            plt.scatter(control[:, 0], control[:, 1])
        if show:
            plt.show()


class Interpolation:
    def __init__(self, knot, inter_points, k=3):
        self.k = k
        self.knot = knot
        self.check_support()
        if type(inter_points) != np.ndarray:
            raise TypeError("We want a np.array")
        if inter_points.shape[1] not in [2, 3]:
            raise ValueError("Wrong shape on ther inter_points")
        self.inter_points = inter_points
        self.N = Spline(self.knot, k=k)
        self.t = self.get_t_points()
        if self.inter_points.shape[0] != len(self.N.N):
            raise Exception("Now we have {} basis functions,and {} interpolation points".format(len(self.N.N),
                            len(inter_points)))
        self.coeff = self.getcoefficients()

    def check_support(self):
        def fullSupport(knot, t, deg): return knot[deg] <= t <= knot[-deg-1]
        def overSupport(knot, t, deg): return (knot[0] == knot[deg+1] or knot[-1] == knot[-deg-2])
        for i in range(len(self.knot)-1):
            if not fullSupport(self.knot, self.knot[i], self.k):
                raise(Exception("You do not have full support on point " + str(self.knot[i]) +" index "+ str(i)))
            elif overSupport(self.knot, self.knot[i], self.k):
                raise(Exception("You have over support dont be stupid"))
        return True

    def get_t_points(self):
        a, b, s = self.knot[0], self.knot[-1], len(self.N.N)-1
        ret = list(a + k * (b-a) / s for k in range(len(self.N.N)))
        return ret

    def __get_a_matrix(self):
        """creats a square tridiagonal matrix of size n"""
        n = len(self.N.N)
        A = np.zeros((n, n))
        t = self.t

        for col in range(n):
            for row in range(n):
                A[row, col] = self.N.eval_basis(t[row], col)
        return A

    def getcoefficients(self):
        self.A = self.__get_a_matrix()
        xcoeff = self.__solve(self.A,self.inter_points[:,0]).reshape((len(self.N.N)))
        ycoeff = self.__solve(self.A,self.inter_points[:,1]).reshape((len(self.N.N)))
        xy = np.empty_like(self.inter_points)
        if self.inter_points.shape[1] == 3:
            zcoeff = self.__solve(self.A,self.inter_points[:,1]).reshape((len(self.N.N)))
            xy[:,2] = zcoeff

        xy[:, 1] = ycoeff
        xy[:, 0] = xcoeff
        return(xy)
    def __solve(self,a,b):
        return np.linalg.solve(a,b.T)

    def x(self, x):
        return (self.splinex.evalfull(x))

    def y(self, y):
        return (self.spliney.evalfull(y))

    def plotbasis(self, show=False):
        self.N.basisplot()
        if show:
            plt.show()

    def plotinter(self, show=False, pts=100):
        self.N.plot2D(self.coeff, self.inter_points)
        if show:
            plt.show()


class Surface:
    def __init__(self, controlnet, knotsu, knotsv: np.dtype, deg: list, weights = None):
        self.ctrnet = controlnet
        self.k = deg
        self.knotsu = knotsu
        self.knotsv = knotsv
        self.splu = Spline(knotsu, deg[0])
        self.splv = Spline(knotsv, deg[1])


        self.weights = weights

    def __call__(self, u, v):
        innerctrpts = np.array(list(self.splu.eval_vector(self.ctrnet[i], u)
                                for i in range(self.ctrnet.shape[1])))
        ret = self.spl.eval_vector(innerctrpts, v)
        return(ret)
        
    def nurbescall(self, u, v):
        
        # ugly solution
        den = 0
        for p in range(len(self.splu.N)):
            evalu = self.splu.eval_basis(u, p)
            for q in range(len(self.splv.N)):
                evalv = self.splv.eval_basis(v, q)
#                print(evalu*evalv*self.weights[q])
                den += evalu*evalv*self.weights[p,q]
        if den == 0:
            raise ValueError("")
        def nurbbasis(i,j,u,v,denominator): return (self.splu.eval_basis(u,i) * 
                      self.splv.eval_basis(v,j) * self.weights[i,j]/denominator)
        ret = np.zeros((3,))
        for i in range(len(self.splu.N)):
            for j in range(len(self.splv.N)):
#                print(nurbbasis(i,j,u,v,den),self.ctrnet[i,j])
#                print(self.ctrnet[i,j].shape)
                ret += nurbbasis(i,j,u,v,den)*self.ctrnet[i,j]
        return ret

    def plotspline(self, pts=200,color="black"):
        if self.weigts == None:
            raise ValueError("You probs wanna do a nurbsplot")
        X, Y = np.meshgrid(np.linspace(self.knotsu[0], self.knotsu[-1], pts),
                            np.linspace(self.knotsv[0], self.knotsv[-1], pts))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for col in range(pts):

            tmp = []
            tmp2 = []
            for row in range(pts):
                ret = self.__call__(X[row, col], Y[row, col])
                ret2 = self.__call__(X[col,row], Y[col,row])
                tmp.append(ret)
                tmp2.append(ret2)
            tmp2 = np.array(tmp2)
            tmp = np.array(tmp)

            ax.plot(tmp[:, 0], tmp[:, 1], tmp[:, 2], alpha=0.5, c=color)
            ax.plot(tmp2[:, 0], tmp2[:, 1], tmp2[:, 2], alpha=0.5, c=color)
    def plotnurbs(self, pts= 200, color="black"):
#        print(len(self.splu.N),len(self.weights), self.ctrnet.shape)
        X, Y = np.meshgrid(np.linspace(self.knotsu[0], self.knotsu[-1], pts),
                            np.linspace(self.knotsv[0], self.knotsv[-1], pts))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for col in range(pts):
            print(col/pts)
            tmp = []
            tmp2 = []
            for row in range(pts):
                ret = self.nurbescall(X[row, col], Y[row, col])
                ret2 = self.nurbescall(X[col,row], Y[col,row])
                tmp.append(ret)
                tmp2.append(ret2)
            tmp2 = np.array(tmp2)
            tmp = np.array(tmp)
#            print(tmp,"tmp \n")
            ax.plot(tmp[:, 0], tmp[:, 1], tmp[:, 2], alpha=0.5, c=color)
            ax.plot(tmp2[:, 0], tmp2[:, 1], tmp2[:, 2], alpha=0.5, c=color)
        


def square3d():
    controlnet = np.array([[[-2,-2,1],[1/2, 1/2, 1],[2,2,1]],
                            [[-2,-2,1], [1, 0, 0],[2,2,1]],
                            [[-2,-2,1],[0, 1, 0],[2,2,1]],
                            [[-2,-2,1],[1/2, 1/2, 1],[2,2,1]]])
    knotsu = knotsv = [0,0,0.5,1,1]
    sur2 = Surface(controlnet,knotsv,knotsu,1)
    sur2.plot()
    plt.savefig("task4.pdf")
    plt.show()




def cone():
    knots   = np.array([0, 0, 0, 0,
                        1/6, 2/6, 3/6, 4/6,
                        1, 1, 1, 1])

    points  = np.array([ [0, 0, 0],
                         [0, 0, 1],
                         [0, 0, 2],
                         [0, 0, 3],
                         [0, 0, 4],
                         [0, 0, 5],
                         [0, 0, 6],
                         [0, 0, 7],
                         [0,0,8]
                          ])

    def func(ctr, r, hm):
        hm1 = hm -1
        A = np.zeros( (hm,hm,3) )
        def a(j): return (j/hm1)*2*np.pi
        for i in range(hm):
            w = 1 - i/hm1
            for j in range(hm):
                A[i, j, 0] = r*(hm1-ctr[i,2])/hm1*np.cos(a(j))
                A[i, j, 1] = r*(hm1-ctr[i,2])/hm1*np.sin(a(j))
                A[i, j, 2] = ctr[i, 2]

        return A
    knots   =[0, 0, 0, 0, 1/6, 2/6, 3/6, 4/6, 5/6 , 1, 1, 1, 1]
    controlnet = func(points, 0.1,9)
    for i in range(controlnet.shape[0]):
        S = Spline(knots, k=3, coeff=controlnet[i])
        I = Interpolation(knots, controlnet[i,:,0:2], k=3)
    #    I.plotinter()
    #    S.plot2D()
        controlnet[i,0:9,0:2] = I.coeff[:,:]

    sur = Surface(controlnet, knots, knots, deg=3)
    sur.plot(50)
    plt.show()

def circle():
    control = np.array([[1, 0], [1, 1],
                        [0, 1], [-1, 1],
                        [-1, 0], [-1, -1],
                        [0, -1], [1, -1],
                        [1, 0]])
    knots = [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]
    weights = np.array([1, np.sqrt(2)/2,
                          1, np.sqrt(2)/2,
                          1, np.sqrt(2)/2,
                          1, np.sqrt(2)/2,
                          1])

    s = Spline(knots, k=[3,3])
    s.plotnurbs(control, weights)
    plt.axis("equal")
    plt.show()
    print(weights)


#def cone2():
control = np.array([[1, 0], 
                    [1, 1],
                    [0, 1], 
                    [-1, 1],
                    [-1, 0], 
                    [-1, -1],
                    [0, -1], 
                    [1, -1],
                    [1, 0]])

weights = np.array([[1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1],
                    [1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1],
                    [1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1],
                    [1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1],
                    [1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1],
                    [1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1],
                    [1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1],
                    [1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1],
                    [1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1, np.sqrt(2)/2, 1]])
A = list()
for j in (range(control.shape[0])):
    A1 = []
    for i in range(control.shape[0]):
        print(j)
        a = np.zeros((3,))
        a[0:2] = (control.shape[0]-1-j)/(control.shape[0]-1) * control[i]
        a[2] = j / (control.shape[0]-1)

        A1.append(a)

    A.append(A1)
    
A = np.array(A)

#controlnet = np.array(list(i/9*control for i in range(9)))

knots = [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]

sur = Surface(A,knots,knots,deg = [2,2],weights = weights)

sur.plotnurbs(pts=100)

plt.show()


    
    
    
    
