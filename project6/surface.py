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

class BilinearPatch:
    def __init__(self, points, domain=(0,1)):
        pass


p1 = array([ [0,0,0], [1,0,1], [2,0,1], [3,0,0] ])
domain = (0, 1)

C = BilinearPatch(p1, domain)


# plt.show()
