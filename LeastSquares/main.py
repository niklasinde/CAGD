from numpy import *
from numpy.linalg import *
from scipy.special import binom
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

class PolynomialEvaluation:
    def __init__(self, points, domain=(0, 1)):
            self.points, self.n, self.domain = points, len(points)-1, domain

    def Casteljau(self, specificBasis=False):
        a, b = self.domain
        c1   = lambda t: (b - t)/b if b != 0. else 0.
        c2   = lambda t: (t - a)/(b - a) if b - a != 0. else 0.

        def B(i:'iteration', p:'degree') -> 'func':
            if p == 0: return lambda _: self.points[i, :]

            return lambda t: c1(t)*B(i, p-1)(t) + c2(t)*B(i+1, p-1)(t)

        return B(0, self.n) if specificBasis == False else B(0, specificBasis)


class Bezier:
    def __init__(self):
        pass

    def __call__(self, i, p):
        return self.getBasis(i, p)

    def getBasis(self, iteration, degree) -> 'func':
        def B(i, p):
            coefficients = binom(p, i)

            def BinomialCoefficients():
                var  = lambda x: (x**i)*((i - x)**(p-i))
                return lambda x: coefficients*var(x)

            return BinomialCoefficients()

        return B(iteration, degree)


class EquationSystem:
    def __init__(self, data, DesiredSize):
        ''' Generates the datavector and basismatrix,
            in the least squares equation system:
                M^t * P = M^t * M * B,
            where B is solved for.
        '''
        self.data, self.n, self.K = data, DesiredSize, len(data)-1

    def __call__(self):
        return self.construct()

    def datavector(self, matrix):
        return dot(matrix, self.data)

    def basismatrix(self):
        B = Bezier()

        C = lambda u, v: list(list( B(iter1, self.n)(u)*B(iter2, self.n)(v)
                                    for iter1 in range(self.n) )
                                    for iter2 in range(self.n)  )

        return lambda U, V: array(list( C(U[k], V[k]) for k in range(self.K) ))

    def construct(self):
        def subdivide(m):
            _m = transpose(m)
            return self.datavector(_m), dot(_m, m)

        U = array([ [0, 0],
                    [1, 1],
                    [1, 0],
                    [0, 1],
                    [0.5, 0.5]
                    ])

        U = linspace(0, 1, self.K)
        return list( subdivide(m) for m in self.basismatrix()(U, U) )


class evalaute:
    def __init__(self, data, DesiredSize):
        Generate = EquationSystem(data, DesiredSize)

        self.b = list( solve(m, p) for m, p in Generate() )

    def surface(self):
        Basis = PolynomialEvaluation(self.b)
        B = Basis.Casteljau()



if __name__=='__main__':
    points = array([ [-2, -2, 1],
                     [2, 2, 1],
                     [1, 0, 0],
                     [0, 1, 0],
                     [1/2, 1/2, 1]
                    ])

    C = evalaute(points, 4)
    C.surface()


