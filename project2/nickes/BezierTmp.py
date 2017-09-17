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
        self.coeff = binom(self.n, self.i)
        self.evaluate = self.__getB()

    def __call__(self, arg:'float') -> 'float':
        return self.evaluate(arg)

    def __getB(self) -> 'float':
        var = lambda x: x**self.i * (1 - x)**(self.n - self.i)
        return lambda arg: self.coeff*var(arg)

class RecurBasis:
    '''Used for task 3, to mimic the
       computation of homework 1.'''
    def __init__(self, pts):
        self.pts, self.n = pts, len(pts)-1
        self.pts2 = self.pts.copy()
        print(self.pts)
        for k in range(self.n + 1):
            tmp = self.domaintransform(self.pts[k, 0])
            print(tmp)
            self.pts[k, 0] = tmp
            print(self.pts[k, 0])
        print(self.pts)
        self.b = self.__getCurve()

    def __getCurve(self) -> 'func':
        def b(i, k):
            if k == 0: return lambda _: self.pts[i, :]
            else: return lambda t: (1-t)*b(i, k-1)(t) + t*b(i+1, k-1)(t)

        return b(0, self.n)

    def domaintransform(self, pt):
        a = (pt - min(self.pts2[:, 0]))
        b = (max(self.pts2[:, 0]) - min(self.pts2[:, 0]))
        return a/b

    def Backwardsdomaintransform(self, pt):
        return ((max(self.pts[:, 0]) - min(self.pts[:, 0]))/
                (pt - min(self.pts[:, 0])))


    def subdivision(self, t=0.5):
        def b(i, k):
            if k == 0: return lambda _: self.pts[i, :]
            else: return lambda t: (1-t)*b(i, k-1)(t) + t*b(i+1, k-1)(t)

        pts1 = array(list(reversed(list(b(0, tmp)(t)
                     for tmp in reversed(range(self.n+1))))))
        pts2 = array(list(b(self.n - tmp, tmp)(t)
                     for tmp in reversed(range(self.n+1))))

        scatter(pts1[:, 0], pts1[:, 1], c='r')
        scatter(pts2[:, 0], pts2[:, 1], c='b')

        A1 = BernPoly(pts1)
        A1.render([-1,t], colour = 'r')
        A2 = BernPoly(pts2)
        A2.render([t, 2], colour = 'b')

    def render(self, nsp=100, colour='r', env=None) -> 'None':
        values = array(list(self.b(x) for x in linspace(0, 1, nsp)))

        plot(values[:, 0], values[:, 1], c=colour)
        if env:
            LineSegment = lambda i, j: [self.pts[i, j], self.pts[i+1, j]]
            scatter(self.pts[:, 0], x[:, 1], c='k', alpha=0.5)
            plot(list(LineSegment[i, 0] for i in range(self.n)),
                 list(LineSegment[i, 1] for i in range(self.n)),
                 c='k', alpha=0.2)

class BernPoly:
    def __init__(self, pts:'list/vector'):
        self.pts, self.n = array(pts), len(pts) - 1
        self.basis = self.__getBasis()
        self.B     = self.__getPolynomial()

    def __call__(self, x) -> 'matrix':
        return self.B(x)

    def __getPolynomial(self) -> 'func':
        return lambda x: sum(self.pts[i]*base(x)
                             for i, base in enumerate(self.basis))

    def __getBasis(self) -> 'matrix':
        return array(list(BernBasis(self.n, idx) for idx in range(self.n + 1)))

    def render(self, domain=[0,1], nsp=100, colour='r', env=None) -> 'None':
        sample = linspace(domain[0], domain[1], nsp)
        u = lambda a: a
        values = array(list(self.B(u(x))
                       for x in linspace(0, 1, nsp)))

        plot(sample, values[:, 1], c=colour, alpha=0.5)

        if env:
            LineSegment = lambda i, j: [self.pts[i, j], self.pts[i+1, j]]
            scatter(self.pts[:, 0], self.pts[:, 1], c='k', alpha=0.5)
            plot(list(LineSegment(i, 0) for i in range(self.n)),
                 list(LineSegment(i, 1) for i in range(self.n)),
                 c='k', alpha=0.2)

    def renderBasisFuncs(self, pt = False) -> 'None':
        tmp = self.basis.copy()
        sample = linspace(0, 1, 100)
        for i, p in enumerate(tmp):
            plot(sample, vectorize(p)(sample))

        xlim(-0.2, 1.2)
        ylim(-0.2, 1.2)

        if pt: print(sum(base(pt) for base in tmp))


if __name__=='__main__':
    x = array([[0,0],
               [0.2, 0.5],
               [0.4, 0.2],
               [0.7, -0.2],
               [1, 0.1]])
    '''task1'''
    #P = BernPoly(x)
    #P.renderBasisFuncs(0.2)
    '''task2'''
    pts3 = array([[0.05, 0.02], [0.1, 0.2],
                 [0.2, -0.1], [0.3, 0],
                 [0.4, 0.1], [0.7, 0.2]])
    #@timeit
    #def RecursiveMethod():
    #    proj1 = RecurBasis(pts)
    #    proj1.render()
    #@timeit
    #def BernsteinMethod():
    #    proj2 = BernPoly(pts)
    #    proj2.render()

    #RecursiveMethod()
    #BernsteinMethod()
    '''task3'''
    pts = array([[-1, 0], [0, 1], [1, -2], [2, 0]])
    B = BernPoly(pts)
    B.render(domain=[-1,2], env=True)
    A = RecurBasis(pts)
    A.subdivision(t=0.4)

    #xlim(-0.2, 1.4)
    #ylim(-0.2, 0.4)
    grid()
    #savefig('SubDivision.pdf')
    show()
