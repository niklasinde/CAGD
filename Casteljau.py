from matplotlib.pyplot import *
from numpy import *
rcParams["figure.figsize"] = [16, 2]

class casteljau():
    def __init__(self, pts:'list[[x, y], []]', *kwargs):
        '''After constructing the object
           with the control points, execute
           the following methods
           self.run(domain[a,b], nsp:=100)
           self.render()'''
        self.pts = pts

    def _b(self, i, k):
        if k == 0:
            return lambda t: self.pts[i, :]
        else:
            return lambda t: (1-t)*self._b(i, k-1)(t) + t*self._b(i+1, k-1)(t)

    def _verify(self, x):
        b = self.pts
        left, right = min(b[:, 0]), max(b[:, 0])
        lower, upper = min(b[:, 1]), max(b[:, 1])
        if left > x[0] or right < x[0]: return False
        elif lower > x[1] or upper < x[1]: return False
        else: return True

    def render(self, *kwargs) -> 'graph':
        b       = self.pts
        coords  = self.coords
        outside = array([pt for pt in coords if self._verify(pt) == False])
        inside  = array([pt for pt in coords if self._verify(pt) == True])
        try:
            plot(outside[:, 0], outside[:, 1], '--r')
        except Exception:
            pass

        try:
            plot(inside[:, 0], inside[:, 1], c='r')
        except Exception:
            pass

        scatter(b[:, 0], b[:, 1])

        for idx in range(len(self.pts)-1):
            plot([b[idx, 0], b[idx+1, 0]], [b[idx, 1], b[idx+1, 1]], c='k')

    def run(self, domain:'list [a, b]', nsp = 100, *kwargs) -> 'self.coords':
        b      = self.pts
        B      = self._b(0, len(self.pts)-1)
        sample = linspace(domain[0], domain[1], nsp)
        self.coords = array([list(B(_)) for _ in sample])


pts = array([[0.05, 0.02], [0.1, 0.2], [0.2, -0.1], [0.3, 0], [0.4, 0.1], [0.7, 0.2]])
p = casteljau(pts)
p.run([0, 1.2])
p.render()
grid()
#savefig('casteljau1.pdf')
show()
