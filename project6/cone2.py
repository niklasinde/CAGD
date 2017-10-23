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

    def evaluate(self, controlpoints, nsp=20):
        N           = self.getBasisGenerator()
        basisvector = array(list( N(i, self.p)
                                  for i in range(len(self.u) - self.p - 1) ))

        getInnerPoints = lambda x: sum( controlpoints[i]*basis(x)
                                        for i,basis in enumerate(basisvector) )

        getValue       = lambda x, y: sum( getInnerPoints(x)*basis(y)
                                           for i,basis in enumerate(basisvector) )

        sample = linspace(self.u[0], self.u[-1], nsp)
        return array(list( list(getValue(x, y) for x in sample ) for y in sample ))



    def NURBS(self, weights):
        def step1(points, w=weights):
            return  array(list( w[i]*append(x, 1)
                               for i,x in enumerate(points) ))

        def step2(obj, tmp):
            return obj.evaluate(controlpoints = tmp)

        def step3(fx):
            return array(list( y[:-1]*y[-1] for y in fx ))

        for points in self.b:
            extended_points = step1(points)
            tmp_object      = Bspline(extended_points, self.u, self.p)

            Yx = step2(tmp_object, extended_points)

            list( ax.plot(Y[:, 0], Y[:, 1], Y[:, 2], c='k', alpha=0.6)
                  for Y in Yx )

        for i, pts in enumerate(self.b):
            points = self.b[:, i, :]
            w2 = array([ 1, 1, 1, 1,
                         1, 1, 1, 1, 1])

            extended_points = step1(points, w2)
            tmp_object      = Bspline(extended_points, self.u, self.p)

            Xy = step2(tmp_object, extended_points)

            #list( ax.plot(Y[:, 0], Y[:, 1], Y[:, 2], c='k', alpha=0.6)
            #      for Y in Xy )

knots   = np.array([0, 0, 0, 0,
                    1/5, 2/5, 3/5, 4/5,
                    1, 1, 1, 1])

points  = np.array([ [0, 0, 0],
                     [0, 0, 1],
                     [0, 0, 2],
                     [0, 0, 3],
                     [0, 0, 4],
                     [0, 0, 5],
                     [0, 0, 6],
                     [0, 0, 7],
                      ])

weights = array([ 0.8, 1, 1, 1,
                  1, 1, 1, 0.8])

def func(ctr, r):
    A = zeros( (8,8,3) )
    def a(j): return (j/7)*2*np.pi
    for i in range(8):
        w = 1 - i/7
        for j in range(8):
            A[i, j, 0] = w*r*sin(a(j))
            A[i, j, 1] = w*r*cos(a(j))
            A[i, j, 2] = ctr[i, 2]

    return A

knots   = np.array([0, 0, 0,
                     1/4, 1/4,
                     2/4, 2/4,
                     3/4, 3/4,
                    1, 1, 1])


points  = np.array([ [0, 0, 0],
                     [0, 0, 1],
                     [0, 0, 2],
                     [0, 0, 3],
                     [0, 0, 4],
                     [0, 0, 5],
                     [0, 0, 6],
                     [0, 0, 7],
                     [0, 0, 8]
                      ])

weights = array([ 1, 1.2, 1,
                  0.98, 1, 0.98,
                  1, 0.98, 1])

weights = array([ 1, 1, 1,
                  1, 1, 1,
                  1, 1, 1])

def func(ctr, r):
    A = zeros( (9,9,3) )

    points  = np.array([ [0, -0.5, 0],
                         [-0.5, -0.5, 1],
                         [-0.5, 0, 2],
                         [-0.5, 0.5, 3],
                         [0, 0.5, 4],
                         [0.5, 0.5, 5],
                         [0.5, 0, 6],
                         [0.5, -0.5, 7],
                         [0, -0.5, 8]
                          ])


    for i in range(9):
        w = 1 - i/8
        A[i, :, :] = points
        for j in range(9):
            A[i, j, 0] = w*A[i, j, 0]
            A[i, j, 1] = w*A[i, j, 1]
            A[i, j, 2] = ctr[i, 2]
            pass

    return A



#for cluster in func(points, 0.1):
#    for pt in cluster: ax.scatter(*pt, c='k')

C = Bspline(func(points, 1), knots, 2)
#print(func(points,1))
C.NURBS(weights)

#for args in points: ax.scatter(*args, c='r')
plt.show()
