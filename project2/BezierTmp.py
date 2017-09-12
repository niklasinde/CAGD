from numpy import *
from scipy.special import binom
from matplotlib.pyplot import *
from functools import reduce
import time

def timeit(method):
    def timed(*args, **kw):
        ts     = time.time()
        result = method(*args, **kw)
        te     = time.time()

        print('%r (%r, %r) %2.2f sec' % (method.__name__, args, kw, te-ts))
        return result
    return timed


class BernBasis:
    def __init__(self, degree:'int', iteration:'int'):
        self.n, self.i = degree, iteration
        assert self.n >= self.i, 'passed max iter.'

    def __call__(self, arg:'float') -> 'float':
        var = lambda x: x**self.i * (1 - x)**(self.n - self.i)
        return binom(self.n, self.i)*var(arg)

class RecurBasis:
    def __init__(self, pts):
        self.pts, self.n = pts, len(pts)-1
        self.Bezier = self._generateCurve()

    def _generateCurve(self) -> 'func':
        def b(i, k):
            if k == 0: return lambda _: self.pts[i, :]
            else: return lambda t: (1-t)*b(i, k-1)(t) + t*b(i+1, k-1)(t)

        return b(0, self.n)

    def _render(self, nsp=100, colour='r', env=None) -> 'None':
        values = array(list(self.Bezier(x) for x in linspace(0, 1, nsp)))

        plot(values[:, 0], values[:, 1], c=colour)
        if env:
            scatter(self.pts[:, 0], x[:, 1], c='k', alpha=0.5)
            plot(list([self.pts[i, 0], self.pts[i+1, 0]] for i in range(self.n)),
                 list([self.pts[i, 1], self.pts[i+1, 1]] for i in range(self.n)),
                 c='k', alpha=0.2)

class BernPoly:
    def __init__(self, pts:'vector'):
        self.pts, self.n = pts, len(pts) - 1
        self.basis = self._generateBasis()
        self.px    = self._generatePolynomial()

    def __call__(self, x) -> 'matrix':
        return self.px(x)

    def _generatePolynomial(self) -> 'func':
        return lambda x: sum(self.pts[i]*base(x)
                             for i, base in enumerate(self.basis))

    def _render(self, nsp=100, colour='r', env=None) -> 'None':
        values = array(list(self.px(x) for x in linspace(0, 1, nsp)))

        plot(values[:, 0], values[:, 1], c=colour)

        if env:
            scatter(self.pts[:, 0], x[:, 1], c='k', alpha=0.5)
            plot(list([self.pts[i, 0], self.pts[i+1, 0]] for i in range(self.n)),
                 list([self.pts[i, 1], self.pts[i+1, 1]] for i in range(self.n)),
                 c='k', alpha=0.2)

    def _PlotBasisFuncs(self, pt=False) -> 'None':
        tmp = self.basis.copy()
        sample = linspace(0, 1, 100)
        for i, p in enumerate(tmp):
            plot(sample, vectorize(p)(sample))

        xlim(-0.2, 1.2)
        ylim(-0.2, 1.2)

        if pt: print(sum(base(pt) for base in tmp))

    def _generateBasis(self) -> 'matrix':
        return array(list(BernBasis(self.n, idx) for idx in range(self.n + 1)))


x = array([[0,0], [0.2, 0.5], [0.4, 0.2], [0.7, -0.2], [1, 0.1]])

'''task1'''
#P = BernPoly(x)
#P._PlotBasisFuncs(0.2)
'''task2'''
pts = array([ [0.05, 0.02], [0.1, 0.2],
              [0.2, -0.1], [0.3, 0],
              [0.4, 0.1], [0.7, 0.2]])
@timeit
def RecursiveMethod():
    proj1 = RecurBasis(pts)
    proj1._render()
@timeit
def BernsteinMethod():
    proj2 = BernPoly(pts)
    proj2._render()

RecursiveMethod()
BernsteinMethod()

'''task3'''

#xlim(-0.2, 1.4)
#ylim(-0.2, 0.4)
#grid()
#savefig('basisfuncs.pdf')
#show()
