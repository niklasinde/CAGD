from numpy import *
from matplotlib.pyplot import *


class casteljau:
    def __init__(self, pts):
        self.b = pts; self.n = len(pts)-1
        self.domain = (min(pts[:, 0]), max(pts[:, 0]))

    def run(self, points):
        def N(i:'iter', p:'degree') -> 'func':
            if p == 0: return lambda _: points[i, :]
            else: return lambda t: (1-t)*N(i, p-1)(t) + t*N(i+1, p-1)(t)

        return N(0, self.n)

    def NURBS(self, w, points):
        return array([ [w[i]*arg[0], w[i]*arg[1], w[i]] for i, arg in enumerate(points) ])

    def evaluate(self, w):
        tmp_points = self.NURBS(w, self.b)
        Basis      = self.run(tmp_points)

        sample = linspace(self.domain[0], self.domain[-1], 100)

        y = array([Basis(x) for x in sample])
        y = self.NURBS(1/w, y)
        for arg in y: scatter(*arg[:2])
        for arg in self.b: scatter(*arg)
        show()

if __name__=='__main__':
    pts = array([[0, 0], [0.2, 0.5], [0.6, -0.2], [1, 0]])
    C = casteljau(pts)
    C.evaluate(1)
    C.evaluate(0.2)
