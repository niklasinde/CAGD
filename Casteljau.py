from matplotlib.pyplot import *
from numpy import *
from itertools import *
rcParams["figure.figsize"] = [16, 2]

class casteljau:
    def __init__(self, pts:'[[x,y],]'):
        '''The valid methods are
           __call__,
           __add__,
           __div__'''
        self.pts = array(pts)
        self.left, self.right = min(self.pts[:, 0]), max(self.pts[:, 0])
        self.lower, self.upper = min(self.pts[:, 1]), max(self.pts[:, 1])
        self.n = len(pts)-1

    def __add__(self, other:'obj'):
        tmp = array([_ for _ in chain(self.pts, other.pts)])
        return casteljau(tmp)

    def __div__(self, other:'obj'):
        pass

    def __call__(self, domain:'list [a, b]', nsp = 100, colour = 'r'):
        B      = self._b(0, self.n)
        sample = linspace(domain[0], domain[1], nsp)
        self.coords = array([list(B(_)) for _ in sample])
        self._render(colour)


    def _b(self, i, k) -> 'func':
        if k == 0:
            return lambda t: self.pts[i, :]
        else:
            return lambda t: (1-t)*self._b(i, k-1)(t) + t*self._b(i+1, k-1)(t)

    def _verify(self, x) -> 'Bolean':
        if x[0] < self.left or x[0] > self.right: return False
        elif x[1] < self.lower or x[1] > self.upper: return False
        else: return True

    def _render(self, colour) -> 'graph':
        coords  = self.coords
        outside = array([pt for pt in coords if self._verify(pt) == False])
        inside  = array([pt for pt in coords if self._verify(pt) == True])
        try:
            plot(outside[:, 0], outside[:, 1], '--r')
        except Exception:
            pass

        try:
            plot(inside[:, 0], inside[:, 1], c=colour)
        except Exception:
            pass

        scatter(self.pts[:, 0],
                self.pts[:, 1])

        for idx in range(len(self.pts)-1):
            plot([self.pts[idx, 0], self.pts[idx+1, 0]],
                 [self.pts[idx, 1], self.pts[idx+1, 1]], c='k')


pts = [[0.05, 0.02], [0.1, 0.2],
       [0.2, -0.1], [0.3, 0],
       [0.4, 0.1], [0.7, 0.2]]
p = casteljau(pts)
p([0, 1.2])

pts2 = [[0.7,0.2], [0.9, -0.1], [1.3, 0]]
g = casteljau(pts2)
g([0, 1.2], colour='b')

f = p + g
f([0,1], colour='k')

grid()
#savefig('casteljau1.pdf')
show()
