from numpy import*
from numpy.linalg import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import coonspatch

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


class BilinearPatch:
    def __init__(self, data, domain=(0,1)):
        self.data = data

    def tmp(self):
        '''
           (1) Construct 2x2 point matrix,
           (2) Compute S(u,v) = Sum(Sum(Bi(u)) Bj(v) bij) linear interp,
           (2.5) in matrix form.
        '''
        def wrapper(points):
            X, Y, Z, C = points

            def getPoint(u, v):
                l = array([ 1-u, u ])
                M = lambda x,y,z,c: array([ [x, y], [z, c] ])
                r = array([ [1-v], [v] ])

                return list( list(dot(dot(l, M( X[i], Y[i], Z[i], C[i] )), r))
                             for i in range(len(X)) )

            return getPoint

        tmpPoints = array([ [0,0,0], [2,0,0], [0,2,0], [0,0,2] ])
        C = wrapper(tmpPoints)

        nsp=3
        sample = linspace(0, 1, nsp)
        XY = array(list( list(C(x, y) for x in sample) for y in sample))
        for X in XY:
            for arg in X:
                #ax.plot(*arg, c='r')
                print(transpose(arg))

            print('--')

            #ax.plot(*X, c='r')


        #ax.scatter(*tmpPoints, c='r')

        #return tmpPoints





    def render(self):
        for arg in self.data: ax.scatter(*arg, c='k')

if __name__=='__main__':
    p1 = array([ [2,2,1], [2,2,1], [1,0,0], [0,1,0], [1/2, 1/2, 1] ])
    domain = (0, 1)

    C = BilinearPatch(p1, domain)
    C.tmp()
    C.render()

    #plt.show()
