from numpy import *
from matplotlib.pyplot import *

class approximate:
    def __init__(self, points, knots, degree):
        self.b, self.u, self.p = points, knots, degree

    def evaluate(self, old_knots, old_points, t=0.5):
        AddKnot = lambda t:  sort(list(append(old_knots.copy(), t)))
        new_knots = AddKnot(t)

        idxset  = lambda t: where(new_knots == t)[0]

        def getPosition(t) -> 'j, h':
            tmp = idxset(t)
            return tmp[-1], len(tmp)

        j, h = getPosition(t)


        def Cj(j, h, u:'knots', b:'points') -> 'func':
            r  = lambda t: t - u[j]
            l  = lambda _: u[j+p - (h-1)] + u[j]
            l2 = lambda _: u[j+p - (h-1)] - u[j]
            c1 = lambda t: 1 - r(t)/l(t) if l(t) != 0 else 0
            c2 = lambda t: r(t)/l2(t) if l2(t) != 0 else 0

            return lambda t: c1(t)*b[j-1] + c2(t)*b[j]

        newpt = Cj(k, h, new_knots, old_points)
        print(newpt)


if __name__=='__main__':
    pts = array([[0, 0], [0.2, 0.4], [0.4, -0.2], [0.6], [1, 0]])
    knots = array([0, 0, 0, 0.5, 1, 1, 1])
    A = approximate(pts, knots, 2)
    A.evaluate(A.b, A.u)
