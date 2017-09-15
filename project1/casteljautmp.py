from matplotlib.pyplot import *
from numpy import *
from itertools import *
from types import *
import contextlib
from prodjec1 import convexhull

@contextlib.contextmanager
def ignored(*exceptions):
    try:
        yield
    except exceptions:
        pass


class casteljau:
    def __init__(self, pts:'[[x,y],]', Basis=None):
        '''Valid methods are: __mul__, __add__, __call__'''
        self.pts, self.n = array(pts), len(pts) - 1
        self.hull = convexhull(list(list(x) for x in self.pts))
        self.__generateBezier(Basis)


    def __mul__(self, other:'obj', Basis=None) -> 'obj':
        return casteljau(array([pt for pt in chain(self.pts, other.pts)]),
                         Basis)


    def __add__(self, other:'obj') -> 'obj':
        def CompositeBezier(t):
            curve1 = self.Bezier(t); curve2 = other.Bezier(t - 1)

            if self.hull(list(curve1)): return curve1
            elif other.hull(list(curve2)): return curve2
            else: return array([0, 0])

        return self.__mul__(other, CompositeBezier)


    def __call__(self, domain:'[l,r]', nsp = 100, colour = 'r') -> 'None':
        sample = linspace(domain[0], domain[1], nsp)
        self.coords = iter(map(self.Bezier, sample))
        self.render(colour)


    def __generateBezier(self, Basis) -> 'func':
        def b(i, k):
            if k == 0: return lambda _: self.pts[i, :]
            else: return lambda t: (1-t)*b(i, k-1)(t) + t*b(i+1, k-1)(t)

        if Basis == None:
            self.Bezier = b(0, self.n)
        else:
            self.Bezier = Basis

    def render(self, colour='r') -> 'None':
        '''Constructs a matplotlib.pyplot object'''
        def PlotCurve(_outside, _inside):
            with ignored(Exception):
                plot(_outside[:, 0], _outside[:, 1], '--k')
            with ignored(Exception):
                plot(_inside[:, 0], _inside[:, 1], c=colour)

        outside, inside = [], []
        for tmp in self.coords:
            '''As of now, convexhull.hull
               requires a list as input.'''
            if self.hull(list(tmp)):
                inside.append(tmp)
            else:
                outside.append(tmp)
        else:
            PlotCurve(array(outside), array(inside))
            self.hull.plot()


if __name__=='__main__':
    pts = array([[0, 0], [0.2, 0.5],
                 [0.4, 0.2], [0.7, -0.2],
                 [1, 0.1]])

    p = casteljau(pts)
    p([0, 1], colour = 'k')

    grid()
    ylim(pts[0, 1] - 0.2, pts[-1, 1] + 0.5)
    xlim(pts[0, 0] - 0.2, pts[-1, 0] + 0.5)
    show()
