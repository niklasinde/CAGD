from numpy import *
from matplotlib.pyplot import *

class Bezier:
    def __init__(self, pts):
        self.pts, self.n = pts, len(pts) - 1
        self.b = self.getB()

    def getB(self):
        c1 = lambda t: (1 - t)
        c2 = lambda t: t

        def b(i, k):
            if k == 0: return lambda _: self.pts[i, :]
            else: return lambda t: c1(t)*b(i, k-1)(t) + c2(t)*b(i+1, k-1)(t)

        return b(0, self.n)

    def render(self):
        values = array((list(self.b(x)
                             for x in linspace(0,
                                               1, 100))))
        plot(values[:, 0], values[:, 1], '--r')
        scatter(self.pts[:, 0], self.pts[:, 1], c= 'k')

class Bspline:
    def __init__(self, pts, knots, deg=2):
        self.pts, self.knots, self.m, self.deg = pts, knots, len(knots) - 1, deg
        if len(self.pts) + self.deg != self.m:
            print('woops', self.m)
            #self.m = self.m - (self.m - (len(self.pts) + self.deg))
            print(self.m)
        self.b = lambda i: self.getB(i)

    def getB(self, idx):
        x = self.knots.copy()

        def N(i, p):
            if p == 0:
                return lambda u: 1 if x[i] <= u <= x[i+1] else 0

            else:
                a1 = lambda u: u - x[i]; b1 = lambda _: x[i+p] - x[i]
                a2 = lambda u: x[i+p+1] - u; b2 = lambda _: x[i+p+1] - x[i+1]

                c1 = lambda u: 0 if abs(b1(u)) < 0.0001 else a1(u)/b1(u)
                c2 = lambda u: 0 if abs(b2(u)) < 0.0001 else a2(u)/b2(u)

                return lambda u: c1(u)*N(i, p-1)(u) + c2(u)*N(i+1, p-1)(u)

        return N(idx, self.deg)

    def render(self):
        C = lambda u: sum(self.pts[i]*self.b(i)(u)
                          for i in range(len(self.pts)))#self.m - self.deg))

        sample = linspace(0, 1, 100)

        X = array([list(C(x)) for x in sample])
        plot(X[:, 0], X[:, 1], '--b')
        scatter(self.pts[:, 0], self.pts[:, 1], c='r')




p = array([[0,0], [0.2, 0.4], [0.5, 1], [0.8, 0.6], [1, 0.3]])
C1 = Bezier(p)
C1.render()

k = array([0, 0, 0.2, 0.2, 0.5, 0.8, 0.85, 0.9, 1.2])

#p2 = array([[0, 0], [0.5, 0.5], [1, 0]])
#print(len(p), len(k), 'n+1', 'p')
C2 = Bspline(p, k, 2)
C2.render()
grid()
show()
