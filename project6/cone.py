from numpy import *
from matplotlib.pyplot import *
from numpy.linalg import solve, norm
#rcParams['figure.figsize'] = 16, 4

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

class Bspline:
    def __init__(self, points, knots, degree):
        self.b, self.u, self.p, self.n = points, knots, degree, len(points)

    def getBasisGenerator(self) -> 'func':
        u = self.u.copy()

        def div(lhs, rhs):
            if rhs == 0: return 0.

            return lhs / rhs

        def recursion(i, k, t):
            if k == 0:
                if u[i] == u[i+1]: return 0.
                if not (u[i] <= t <= u[i+1]): return 0.

                return 1.

            return div( (t - u[i]), (u[i+k] - u[i]) )*recursion(i, k-1, t)\
                   + div( (u[i+k+1] - t), (u[i+k+1] - u[i+1]) )*recursion(i+1, k-1, t)

        def N(i, k) -> 'func':
            return lambda t: recursion(i, k, t)

        return N

    def evaluate(self, controlpoints, nsp=10):
        b           = controlpoints
        N           = self.getBasisGenerator()
        basisvector = array(list( N(i, self.p) for i in range(len(self.u) - self.p - 1) ))

        n = self.n
        a  = list( i/n*2*pi for i in range(n+1) )
        polar = array(list( [sin(alpha), 0, 0, 0] for alpha in a))
        print(polar)

        getInnerPoints = lambda x: sum( (self.b[i]+polar[i])*basis(x)
                                  for i,basis in enumerate(basisvector) )

        getValue = lambda x, y: sum( getInnerPoints(x)*self.b[i]*basis(y)
                                  for i,basis in enumerate(basisvector) )

        sample = linspace(self.u[0], self.u[-1], nsp)
        return array(list( list(getValue(x, y) for x in sample ) for y in sample )),\
               array(list( list(getValue(y, x) for x in sample ) for y in sample ))


    def NURBS(self, weights):
        def step1():
            return  array(list( weights[i]*append(x, 1)
                               for i,x in enumerate(self.b.copy()) ))

        def step2(obj, tmp):
            return obj.evaluate(controlpoints = tmp)

        def step3(fx):
            return array(list( y[:-1]*y[-1] for y in fx ))

        extended_points = step1()
        tmp_object      = Bspline(extended_points, self.u, self.p)

        Horizontal, Vertical = step2(tmp_object, extended_points)
        YX = list( step3( _ ) for _ in Horizontal )
        XY = list( step3( _ ) for _ in Vertical )

        list( ax.plot(Y[:, 0], Y[:, 1], Y[:, 2], c='k', alpha=0.5)
              for Y in YX )
        list( ax.plot(Y[:, 0], Y[:, 1], Y[:, 2], c='k', alpha=0.5)
              for Y in XY)



knots   = np.array([0, 0, 0, 1/2, 1, 1, 1])

points  = np.array([ [0, 0, 0],
                     [0, 0.5, 0],
                     [0, 0.7, 0],
                     [0, 1, 0],
                     ])

weights = array([ 1, 1, 1,
                  1])

C = Bspline(points, knots, 2)

C.NURBS(weights)

for args in points: ax.scatter(*args, c='r')
plt.show()
