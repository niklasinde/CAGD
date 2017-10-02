from numpy import *
from matplotlib.pyplot import *

class approximate:
    def __init__(self, points, knots, degree):
        self.b, self.u, self.p = points, knots, degree

    def evaluate(self):
        AddKnot = lambda t:  sort(append(self.u.copy(), t))
        idx     = lambda t: where(AddKnot(t) == t)[0]
        def analyze(t):
            tmp = idx(t)
            return tmp[-1], len(tmp)

        j, a = analyze(t)
        def Cj(i, h, x, v):
            a = t-
        inner = lambda x, v:  append(x, Cj(i, h, x, v))
        A = iter([inner(*args) for range(a, p + 2)])

        pts = self.b.copy()

        for c in A:
            pts = c(pts)
        else:
            return pts


if __name__=='__main__':
    pts = array([[0, 0], [0.2, 0.4], [0.4, -0.2], [0.6], [1, 0]])
    knots = array([0, 0, 0, 0.5, 1, 1, 1])
    A = approximate(pts, knots, 2)
    A.evaluate()
