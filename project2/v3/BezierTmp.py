from numpy import *
from scipy.special import binom
from matplotlib.pyplot import *
from functools import reduce
from TrivialReject import rectangle
import time
import contextlib

@contextlib.contextmanager
def ignored(*exceptions):
    try:
        yield
    except exceptions:
        pass

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
    def __init__(self, pts):
        self.pts, self.n = pts, len(pts)-1
        self.b       = self.__getCurve(0, self.n)
        self.domain  = [min(self.pts[:, 0]), max(self.pts[:, 0])]
        self.ydomain = [min(self.pts[:, 1]), max(self.pts[:, 1])]

    def __getCurve(self, s, v, env=False) -> 'func':
        if env == False:
            def c1(t, _):
                a = self.domain[1] - t; b = self.domain[1] - self.domain[0]
                return a/b if abs(b) > 0.01 else 0

            def c2(t, _):
                a = t - self.domain[0]; b = self.domain[1] - self.domain[0]
                return a/b if abs(b) > 0.01 else 0
        else:
            def c1(t, i):
                return (i+1)/(self.n+1)
            def c2(t, i):
                return 1 - (i+1)/(self.n+1)

        def b(i, k):
            if k == 0: return lambda _: self.pts[i, :]
            else: return lambda t: c1(t,i)*b(i, k-1)(t) + c2(t,i)*b(i+1, k-1)(t)

        return b(s, v)

    def intersections(self, other:'obj'):
        self.other1 = other
        self.other = rectangle(self.other1.domain, self.other1.ydomain)

        def getSubDivision(obj):
            return obj.SubDivision(t = abs(obj.domain[0] - obj.domain[1])/2)

        def performTrivialReject(A, B):
            a = rectangle(A.domain, A.ydomain)
            b = rectangle(B.domain, B.ydomain)
            if a.control(self.other):
                if b.control(self.other): return A, B
                else: return A
            elif a.control(self.other): return B
            else: raise Exception

        def wrapper(X):
            out = performTrivialReject(*getSubDivision(X))
            return iter(out) if shape(out) != () else [out]

        def Loop(X):
            for A in wrapper(X):
                try:
                    Loop(A)
                except Exception:
                    break

        with ignored(Exception):
            Loop(self)

    def SubDivision(self, t=0.5) -> 'obj, obj':
        pts1 = array(list(reversed(list(self.__getCurve(0, tmp)(t)
                     for tmp in reversed(range(self.n+1))))))
        pts2 = array(list(self.__getCurve(self.n - tmp, tmp)(t)
                     for tmp in reversed(range(self.n+1))))

        return RecurBasis(pts1), RecurBasis(pts2)

    def DegreeElevation(self) -> 'None':
        tmp = [list(self.pts[0].copy())]
        [tmp.append(x)
         for x in list(list(self.__getCurve(i, 1, env=True)(0.5))
                       for i in range(self.n))]

        tmp.append(list(self.pts[-1].copy()))

        tmp = array(tmp)

        self.pts, self.n = tmp, len(tmp)-1
        self.b       = self.__getCurve(0, self.n)
        self.domain  = [min(self.pts[:, 0]), max(self.pts[:, 0])]
        self.ydomain = [min(self.pts[:, 1]), max(self.pts[:, 1])]

    def render(self, nsp=100, colour='r', env=True) -> 'None':
        domain = self.domain
        values = array(list(self.b(x)
                            for x in linspace(domain[0], domain[1], nsp)))

        plot(values[:, 0], values[:, 1], c=colour)

        if env:
            print('env active')
            LineSegment = lambda i, j: [self.pts[i, j], self.pts[i+1, j]]
            scatter(self.pts[:, 0], self.pts[:, 1], c='k', alpha=0.5)
            plot(list(LineSegment(i, 0) for i in range(self.n)),
                 list(LineSegment(i, 1) for i in range(self.n)),
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
        values = array(list(self.B(u(x))
                       for x in linspace(0, 1, nsp)))

        plot(values[: 0], values[:, 1], alpha=0.5, c=colour)

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
    #x = array([[0,0],
    #           [0.2, 0.5],
    #           [0.4, 0.2],
    #           [0.7, -0.2],
    #           [1, 0.1]])
    '''task1'''
    #P = BernPoly(x)
    #P.renderBasisFuncs(0.2)
    '''task2'''
    #pts3 = array([[0.05, 0.02], [0.1, 0.2],
    #             [0.2, -0.1], [0.3, 0],
    #             [0.4, 0.1], [0.7, 0.2]])
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
    #A.SubDivision(t=0.6)
    '''task4'''
    #pts = array([[0,0], [9, -4], [7, 5], [2, -4]])
    #linepts = array([[4,5], [6, -4]])
    #pts = array([[0,0], [0.25, 0.75], [0.5, 1], [0.75, 0.75], [1, 0]])
    #linepts = array([[-0.1, 0.51], [1.1, 0.5]])
    #A = RecurBasis(pts)
    #L = RecurBasis(linepts)
    #A.render(colour = 'r', env=False)
    #L.render(colour = 'r', env=False)
    #A.intersections(L)

    '''task6'''
    #pts = array([[0, -0.1], [0.2, 0.4], [0.5, 0.5], [0.8, 0.3], [0.6, -0.2]])
    pts = array([[0,0], [1,1], [2,1]])
    A = RecurBasis(pts)
    A.render()
    #for _ in range(10):
    #    if _ == 0 or _ == 5: tmp = True
    #    else: tmp = False
    #    A.DegreeElevation()
    #    A.render(env=tmp,colour='k')

    A.DegreeElevation()
    A.render(colour='k')

    pts2 = A.pts.copy()

    pts2[2, :] = pts[1, :]
    pts2[1, : ] = [0.2, 0]
    scatter(pts2[:, 0], pts2[:,1])
    print(pts2)
    B = RecurBasis(pts2)
    B.render(env=False)

    #xlim(-0.2, 1.4)
    #ylim(-0.2, 1.4)
    grid()
    #savefig('DegeeElevation3.pdf')
    show()
