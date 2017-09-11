from matplotlib.pyplot import *
from numpy import *
from itertools import *
from types import *
import contextlib
from prodjec1 import convexhull
rcParams["figure.figsize"] = [16, 2]
@contextlib.contextmanager
def ignored(*exceptions):
    try:
        yield
    except exceptions:
        pass

class casteljau:
    def __init__(self, pts:'[[x,y],]', Basis=None):
        '''The valid methods are
           __mul__,
           __div__,
           __add__,
           __call__'''
        self.pts, self.n = array(pts), len(pts)-1
        self._generateBezier(Basis)
        self._generateHull()

    def __mul__(self, other:'obj', Basis=None):
        assert self._TrivialReject(other)
        return casteljau(array([pt for pt in chain(self.pts, other.pts)]),
                         Basis)

    def __truediv__(self, m:'integer'):
        assert type(m) is int, 'TypeError'
        pass

    def __add__(self, other:'obj'):
        assert self._TrivialReject(other)
        def CompositeBezier(t):
            curve1 = self.Bezier(t)
            curve2 = other.Bezier(t-1)
            if self._assertHull(curve1):
                return curve1
            elif other._assertHull(curve2):
                return curve2
            else:
                return array([0, 0])

        return self.__mul__(other, CompositeBezier)

    def __call__(self, domain:'list [a, b]', nsp = 100, colour = 'r'):
        self.coords = iter(map(self.Bezier
                               ,linspace(domain[0], domain[1], nsp)))
        self._render(colour)

    def _generateBezier(self, Basis) -> 'func':
        def b(i, k):
            if k == 0: return lambda _: self.pts[i, :]
            else: return lambda t: (1-t)*b(i, k-1)(t) + t*b(i+1, k-1)(t)

        if Basis == None:
            self.Bezier = b(0, self.n)
        else:
            self.Bezier = Basis

    def _generateHull(self):
        hull = convexhull(self.pts)
        #self.left, self.right = min(self.pts[:, 0]), max(self.pts[:, 0])
        #self.lower, self.upper = min(self.pts[:, 1]), max(self.pts[:, 1])

    def _assertHull(self, x) -> 'Bolean':
        return hull.checker(x)

        #if x[0] < self.left or x[0] > self.right: return False
        #elif x[1] < self.lower or x[1] > self.upper: return False
        #else: return True

    def _TrivialReject(self, other:'obj'):
        return True

    def _render(self, colour='r') -> 'graph':
        def PlotCurve(_outside, _inside):
            with ignored(Exception):
                plot(_outside[:, 0], _outside[:, 1], '--k')
            with ignored(Exception):
                plot(_inside[:, 0], _inside[:, 1], c=colour)

        def PlotPoints():
            scatter(self.pts[:, 0],
                    self.pts[:, 1], c='k', alpha=0.5)

            for idx in range(self.n):
                plot([self.pts[idx, 0], self.pts[idx+1, 0]],
                     [self.pts[idx, 1], self.pts[idx+1, 1]], c='k', alpha=0.2)

        outside, inside = [], []
        for tmp in self.coords:
            if self._assertHull(array(tmp)):
                inside.append(tmp)
            else:
                outside.append(tmp)
        else:
            PlotCurve(array(outside), array(inside))
            PlotPoints()


pts = [[0.05, 0.02], [0.1, 0.2],
       [0.2, -0.1], [0.3, 0],
       [0.4, 0.1], [0.7, 0.2]]
p = casteljau(pts)
p([0, 1], colour = 'k')


#pts2 = [[0.7,0.2], [0.9, -0.1], [1.3, 0]]
#g = casteljau(pts2)
#g([0, 1], colour='k')


#f = p * g
#f([0,1], colour='k')

#s = p + g
#s([0, 2])

#grid()
#ylim(-0.2, 0.4)
#xlim(-0.2, 1.4)
#savefig('casteljau3.pdf')
#show()
