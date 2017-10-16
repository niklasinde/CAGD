from numpy import*
from numpy.linalg import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

class Bezier:
    def __init__(self, points, domain=(0,1)):
        self.b, self.n = points, len(points)-1
        self.domain = domain
        self.basis  = self.__run()

    def __call__(self, s):
        return self.basis(s)

    def __run(self):
        c1 = lambda t: (self.domain[-1] - t)/self.domain[-1] if self.domain[-1] != 0 else 0
        c2 = lambda t: (t - self.domain[0])/(self.domain[-1] - self.domain[0]) if (self.domain[-1] - self.domain[0]) != 0 else 0

        def B(i, p):
            if p == 0: return lambda _: self.b[i, :]
            return lambda t: c1(t)*B(i, p-1)(t) + c2(t)*B(i+1, p-1)(t)

        return B(0, self.n)

    def evaluate(self, nsp=100):
        basis  = self.basis
        sample = linspace(self.domain[0], self.domain[-1], nsp)

        return array(list( basis(x) for x in sample ))

    def render(self):
        Y = self.evaluate()
        ax.plot(xs = Y[:, 0], ys = Y[:, 1], zs = Y[:,2])

        for arg in self.b: ax.scatter(*arg, c='k')
        ax.plot(self.b[:, 0], self.b[:, 1], self.b[:,2], alpha=0.2, c='k')


class surface:
    def __init__(self, points1, points2, points3, points4, domain=(0,1)):
        ''' points1, points2 are on the same axis, and
            points3, points4 are on the same axis.
        '''
        self.b1, self.b2, self.b3, self.b4 = points1, points2, points3, points4
        self.domain = domain
        _ = self.domain
        self.c0, self.c1, self.d0, self.d1 = Bezier(points1, domain=_),\
                                             Bezier(points2, domain=_),\
                                             Bezier(points3, domain=_),\
                                             Bezier(points4, domain=_)
        self.C = self.__run()

    def __run(self):
        _p = lambda t: (self.domain[-1] - t)/self.domain[-1] if self.domain[-1] != 0 else 0
        _q = lambda t: (t - self.domain[0])/(self.domain[-1] - self.domain[0])\
                       if (self.domain[-1] - self.domain[0]) != 0 else 0


        Lc = lambda s, t: _p(t)*self.c1(s) + _q(t)*self.c0(s)
        Ld = lambda s, t: _p(s)*self.d0(t) + _q(s)*self.d1(t)
        B  = lambda s, t: self.c1(0)*_p(s)*_p(t)\
                          + self.c1(1)*_q(s)*_p(t)\
                          + self.c0(0)*_p(s)*_q(s)\
                          + self.c0(1)*_q(s)*_q(s)

        return lambda s, t: Lc(s, t) + Ld(s, t) - B(s, t)

    def evaluate(self, nsp=30):
        X = linspace(self.domain[0], self.domain[1], nsp)
        Y = linspace(self.domain[0], self.domain[1], nsp)

        XY = lambda _: array(list( list( self.C(x, y) for x in X ) for y in Y ))
        YX = lambda _: array(list( list( self.C(y, x) for x in X ) for y in Y ))

        return XY(0), YX(0)

    def render(self):
        XY, YX = self.evaluate()
        for Y in XY:
            ax.plot(xs = Y[:, 0], ys = Y[:, 1], zs = Y[:,2], alpha=0.5, c='k')
        for Y in YX:
            ax.plot(xs = Y[:, 0], ys = Y[:, 1], zs = Y[:,2], alpha=0.5, c='k')



p1 = array([ [0,0,0], [1,0,1], [2,0,1], [3,0,0] ])
p2 = array([ [0,3,0], [1,3,1], [2,3,1], [3,3,0] ])
p3 = array([ [0,0,0], [0,1,1], [0,2,-3], [0,3,0] ])
p4 = array([ [3,0,0], [3,1,1], [3,2,1], [3,3,0] ])

#p1 = array([ [0,0,0], [1/3,0,1/3], [2/3,0,1/3], [3/3,0,0] ])
#p2 = array([ [0,3/3,0], [1/3,3/3,1/3], [2/3,3/3,1/3], [3/3,3/3,0] ])
#p3 = array([ [0,0,0], [0,1/3,1/3], [0,2/3,-3/3], [0,3/3,0] ])
#p4 = array([ [3/3,0,0], [3/3,1/3,1/3], [3/3,2/3,1/3], [3/3,3/3,0] ])


lista = (p1, p2, p3, p4)
domain = (0, 3)

C = surface(*lista, domain)
C.render()




for arg in lista:
    C = Bezier(arg, domain )
    C.render()

plt.show()
