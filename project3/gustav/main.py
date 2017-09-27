from numpy import *
from matplotlib.pyplot import *

from Basis import *

class Curve:
    def __init__(self, pts, knots, deg):
        self.pts = pts; self.knots = knots; self.n = len(knots) - 1
        self.deg = deg
        self.domain = [min(self.pts[:, 0]), max(self.pts[:, 0])]

        B = Bspline(self.pts, self.knots)
        self.b = lambda i, d: B.getBasis(i, d)

    def evaluate(self):
        C = lambda u: sum(self.pts[i]*self.b(i, self.deg)(u)
                          for i in range(self.n - self.deg))

        sample = linspace(self.domain[0], self.domain[1], 1000)
        X = array([list(C(x)) for x in sample])

        plot(X[:, 0], X[:, 1])
        print(self.knots)

        K = array([list(C(x)) for x in self.knots])
        scatter(self.knots, K[:, 1], c='r')
        grid()
        #show()

if __name__=='__main__':
    controlpts = array([[0, 0], [3, 4],
                        [7, 5], [9, 2], [13, 1], [10, 1], [7,1]])
    #controlpts = array([[0,0], [0.2, 1], [0.4, 0.4], [0.5, 0.3],
    #                    [0.6, 0.25], [0.7, -0.2], [1, 0]])

    #scatter(controlpts[:, 0], controlpts[:, 1])
    knots = array([0,0,0 , 0.5, 0.6, 0.7, 0.8, 1, 1, 1])

    #B = Bspline(controlpts, knots)
    #B.getBasis()
    #A = Curve(controlpts, knots, 2)
    B = Bspline(controlpts, knots)
    B.renderBasis(2)
    #A.evaluate()
    show()

    plsavefig('problem4.pdf')

