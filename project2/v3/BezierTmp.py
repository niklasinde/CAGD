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
        self.b = self.__getCurve(0, self.n)
        self.domain = [min(self.pts[:, 0]), max(self.pts[:, 0])]
        self.ydomain = [min(self.pts[:, 1]), max(self.pts[:, 1])]

    def __getCurve(self, s, v) -> 'func':
        def c1(t):
            a = self.domain[1] - t
            b = self.domain[1] - self.domain[0]
            return a/b if abs(b) > 0.01 else 0

        def c2(t):
            a = t - self.domain[0]
            b = self.domain[1] - self.domain[0]
            return a/b if abs(b) > 0.01 else 0

        def b(i, k):
            if k == 0: return lambda _: self.pts[i, :]
            else: return lambda t: c1(t)*b(i, k-1)(t) + c2(t)*b(i+1, k-1)(t)

        return b(s, v)

    def inside(self, other:'obj'):
        X   = list(self.domain).copy(); [X.append(x) for x in other.domain]
        tmp = X.copy(); tmp.sort()

        Y    = list(self.ydomain).copy(); [Y.append(y) for y in other.ydomain]
        tmpy = Y.copy(); tmpy.sort()

        d = 0.2
        plot(X[0:2], [Y[0], Y[0]], c='k', alpha=d)
        plot(X[0:2], [Y[1], Y[1]], c='k', alpha=d)
        plot([X[0], X[0]], Y[0:2], c='k', alpha=d)
        plot([X[1], X[1]], Y[0:2], c='k', alpha=d)

        if abs(X[0] - X[1]) <= 0.01:
            scatter(self.domain[0], self.ydomain[0], c='r')
            self.render(env=False)
            raise Exception

        def checkBoundary(A, B):
            for a in A:
                if a in set(B): print('in'); return True
            else: print('out'); return False

        def checkInside(A, B):
            C = (not A[0:2] == B[0:2] or not A[0:2] == B[2:4])
            print(C)
            return C

        if checkBoundary(self.domain, other.domain) or checkInside(tmp, X):

            if (checkBoundary(self.ydomain, other.ydomain)
                or checkInside(tmpy, Y)):

                return True

            else: return False
        else: return False


    def intersections(self, other:'obj'):
        def evaluateDivision(obj):
            t = abs(obj.domain[0] - obj.domain[1])/2
            s = obj.b(t)
            scatter(s[0], s[1])
            B1, B2 = obj.subdivision(t)

            if B1.inside(other):
                if B2.inside(other): return B1, B2
                else: return B1
            elif B2.inside(other): return B2
            else: raise Exception

        def wrapper(X):
            tmp = evaluateDivision(X)
            return iter(tmp) if shape(tmp) != () else [tmp]

        def Loop(X):
            for A in wrapper(X):
                try:
                    Loop(A)
                except Exception:
                    #A.render(env=False)
                    #break
                    raise Exception

        self.inside(other) #Temporary
        try:
            Loop(self)
        except Exception:
            pass

    def subdivision(self, t=0.5) -> 'obj, obj':
        pts1 = array(list(reversed(list(self.__getCurve(0, tmp)(t)
                     for tmp in reversed(range(self.n+1))))))
        pts2 = array(list(self.__getCurve(self.n - tmp, tmp)(t)
                     for tmp in reversed(range(self.n+1))))

        A1 = RecurBasis(pts1)
        A2 = RecurBasis(pts2)
        return A1, A2

    def render(self, nsp=100, colour='r', env=True) -> 'None':
        domain = self.domain
        values = array(list(self.b(x)
                            for x in linspace(domain[0], domain[1], nsp)))

        plot(values[:, 0], values[:, 1])#, c=colour)

        if env:
            print('env active')
            LineSegment = lambda i, j: [self.pts[i, j], self.pts[i+1, j]]
            scatter(self.pts[:, 0], self.pts[:, 1], c='k', alpha=0.5)
            plot(list(LineSegment(i, 0) for i in range(self.n)),
                 list(eLineSegment(i, 1) for i in range(self.n)),
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

        plot(sample, values[:, 1], alpha=0.5) #,c=colour

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
    #pts = array([[-1, 0], [0, 1], [1, -2], [2, 0]])
    #B = BernPoly(pts)
    #B.render(domain=[0, 2], env=True)
    #S = RecurBasis(pts)
    #S.render()
    #A = RecurBasis(pts)
    #A.render()
    #A.subdivision(t=0.6)
    '''task4'''
    pts = array([[0,0], [9, -4], [7, 5], [2, -4]])
    linepts = array([[4,5], [6, -4]])
    A = RecurBasis(pts)
    L = RecurBasis(linepts)
    A.render(colour = 'r', env=False)
    L.render(colour = 'r', env=False)
    A.intersections(L)

    #xlim(-0.2, 1.4)
    #ylim(-0.2, 0.4)
    #grid()
    #savefig('SubDivision2.pdf')
    show()
