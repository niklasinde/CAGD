from matplotlib.pyplot import *
from numpy import *
from itertools import *
from types import *
rcParams["figure.figsize"] = [16, 2]

class casteljau:
    def __init__(self, pts:'[[x,y],]', Basis=None):
        '''The valid methods are
           __mul__,
           __div__,
           __add__,
           __call__'''
        self.pts = array(pts)
        self.left, self.right = min(self.pts[:, 0]), max(self.pts[:, 0])
        self.lower, self.upper = min(self.pts[:, 1]), max(self.pts[:, 1])
        self.n = len(pts)-1
        if Basis == None:
            self.B = self._b(0, self.n)
        else:
            self.B = Basis

    def __mul__(self, other:'obj', Basis=None):
        assert self._TrivialReject(other)
        tmp = array([_ for _ in chain(self.pts, other.pts)])
        return casteljau(tmp, Basis)

    def __truediv__(self, m:'integer'):
        assert type(m) is int, 'TypeError'
        pass

    def __add__(self, other:'obj'):
        assert self._TrivialReject(other)
        def test(t):
            curve1 = self.B(t)
            curve2 = other.B(t-1)
            if self._InHull(curve1):
                return curve1
            elif other._InHull(curve2):
                return curve2
            else:
                return array([0, 0])

        tmp = self.__mul__(other, test)
        return tmp

    def __call__(self, domain:'list [a, b]', nsp = 100, colour = 'r'):
        sample = linspace(domain[0], domain[1], nsp)
        self.coords = array([list(self.B(t)) for t in sample])
        self._render(colour)


    def _b(self, i, k) -> 'func':
        if k == 0:
            return lambda t: self.pts[i, :]
        else:
            return lambda t: (1-t)*self._b(i, k-1)(t) + t*self._b(i+1, k-1)(t)

    def _InHull(self, x) -> 'Bolean':
        if x[0] < self.left or x[0] > self.right: return False
        elif x[1] < self.lower or x[1] > self.upper: return False
        else: return True

    def _TrivialReject(self, other:'obj'):
        return True

    def _render(self, colour='r') -> 'graph':
        def PlotCurve(_outside, _inside):
            try:
                plot(_outside[:, 0], _outside[:, 1], '--r')
            except Exception:
                pass

            try:
                plot(_inside[:, 0], _inside[:, 1], c=colour)
            except Exception:
                pass

        def PlotPoints():
            scatter(self.pts[:, 0],
                    self.pts[:, 1], c='k', alpha=0.5)

            for idx in range(self.n):
                plot([self.pts[idx, 0], self.pts[idx+1, 0]],
                     [self.pts[idx, 1], self.pts[idx+1, 1]], c='k', alpha=0.2)

        outside = array([pt for pt in self.coords if self._InHull(pt) == False])
        inside  = array([pt for pt in self.coords if self._InHull(pt) == True])
        PlotCurve(outside, inside)
        PlotPoints()


pts = [[0.05, 0.02], [0.1, 0.2],
       [0.2, -0.1], [0.3, 0],
       [0.4, 0.1], [0.7, 0.2]]
p = casteljau(pts)
#p([0, 1.2])

pts2 = [[0.7,0.2], [0.9, -0.1], [1.3, 0]]
g = casteljau(pts2)
#g([0, 1.2], colour='b')

#f = p * g
#f([0,1], colour='k')

s = p + g
s([0, 2])

grid()
#savefig('casteljau1.pdf')
show()
