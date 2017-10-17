from numpy import*
from scipy.special import binom
from numpy.linalg import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

class Bezier:
    def __init__(self, points, domain=(0,1)):
        self.b, self.n = points, len(points)-1
        self.domain = domain
        self.degree = len(points) - 1

    def __call__(self, s):
        return self.BernPolynomial(s)

    def Bern(self, degree, iteration):
        coeff = binom(degree, iteration)
        def getB():
            var = lambda x: x**iteration * (1-x)**(degree - iteration)
            return lambda arg: coeff*var(arg)

        return getB()

    def BernPolynomial(self, points):
        BASIS = array(list( self.Bern(self.degree, idx)
                            for idx in range(self.degree + 1) ))

        return lambda x: sum( points[i]*basis(x)
                              for i, basis in enumerate(BASIS))

    def evaluate(self, points, nsp=100):
        poly   = self.BernPolynomial(points)
        sample = linspace(self.domain[0], self.domain[-1], nsp)

        return array(list( poly(x) for x in sample ))

    def render(self):
        Y = self.evaluate(self.b)
        ax.plot(xs = Y[:, 0], ys = Y[:, 1], zs = Y[:,2])

        for arg in self.b: ax.scatter(*arg, c='k')
        ax.plot(self.b[:, 0], self.b[:, 1], self.b[:,2], alpha=0.2, c='k')

        XY = self.evaluate(Y)
        ax.plot(XY[:, 0], XY[:, 1], XY[:, 2])

class BilinearPatch:
    def __init__(self, data, domain=(0,1)):
        self.data = data

    def render(self):
        for arg in self.data: ax.scatter(*arg, c='k')

if __name__=='__main__':
    p1 = array([ [2,2,1], [2,2,1], [1,0,0], [0,1,0], [1/2, 1/2, 1] ])
    domain = (0, 1)

    A = Bezier(p1)
    A.render()

    #C = BilinearPatch(p1, domain)
    #C.render()

    plt.show()
