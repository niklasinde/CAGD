from numpy import*
from numpy.linalg import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

class Bezier:
    def __init__(self, points, domain=False):
        self.b, self.n = points, len(points)-1
        self.domain = (min(points[:,0]), max(points[:, 0])) if domain == False else domain

    def __run(self):
        c1 = lambda t: (self.domain[-1] - t)/self.domain[-1] if self.domain[-1] != 0 else 0
        c2 = lambda t: (t - self.domain[0])/(self.domain[-1] - self.domain[0]) if (self.domain[-1] - self.domain[0]) != 0 else 0

        def B(i, p):
            if p == 0: return lambda _: self.b[i, :]
            return lambda t: c1(t)*B(i, p-1)(t) + c2(t)*B(i+1, p-1)(t)

        return B(0, self.n)

    def evaluate(self, nsp=100):
        basis  = self.__run()
        sample = linspace(self.domain[0], self.domain[-1], nsp)

        return array(list( basis(x) for x in sample ))

    def render(self):
        Y = self.evaluate()
        ax.plot(xs = Y[:, 0], ys = Y[:, 1], zs = Y[:,2])

        for arg in self.b: ax.scatter(*arg, c='k')
        ax.plot(self.b[:, 0], self.b[:, 1], self.b[:,2], alpha=0.2, c='k')


class surface:
    def __init__(self, points1, points2, points3, points4):
        ''' points1, points2 are on the same axis, and
            points3, points4 are on the same axis.
        '''
        self.b1, self.b2, self.b3, self.b4 = points1, points2, points3, points4
        self.n = len(self.b1[0, :]); self.m = len(self.b1[:, 0])

    def getPoints(self):
        n, m = self.n-1, self.m-1

        def getbu(M):
            return lambda i, j: (1 - i/m)*M[0, j] + i/m*M[m, j]

        def getbv(M):
            return lambda i, j: (1 - i/m)*M[i, 0] + i/m*M[i, m]

        def getbuv(M):
            def tmp(i,j):
                q = array([1-i/m, i/m])
                w = array([ [M[0,0], M[0,m]], [M[m,0], M[m,m]] ])
                r = array([ [1-i/m], [i/m] ])

                return dot( dot(q, w), r)

            return tmp

        def tmp(axis):
            p1, p2, p3, p4 = self.b1, self.b2, self.b3, self.b4
            M = zeros( (self.m, self.m) )

            M[0,:]  = p1[:, axis]
            M[-1,:] = p2[:, axis]

            M[:,0]  = p3[:, axis]
            M[:,-1] = p4[:, axis]

            bu, bv, buv = getbu(M), getbv(M), getbuv(M)

            def b(i,j):
                M[i, j] = bu(i, j) + bv(i, j) + buv(i, j)

            [[b(j, i) for i in range(1,m)] for j in range(1,m)]

            return M

        points = tmp(idx)
        Bezier(points, domain=(0,3))
        B.render()




p1 = array([ [0,0,0], [1,0,1], [2,0,1], [3,0,0] ])
p2 = array([ [0,3,0], [1,3,1], [2,3,1], [3,3,0] ])
p3 = array([ [0,0,0], [0,1,1], [0,2,-3], [0,3,0] ])
p4 = array([ [3,0,0], [3,1,1], [3,2,1], [3,3,0] ])

lista = (p1, p2, p3, p4)

C = surface(p1, p2, p3, p4)
C.getPoints()

#for arg in lista:
#    C = Bezier(arg, domain=(0, 3) )
#    C.render()

#plt.show()
