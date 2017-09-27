from numpy import *
from matplotlib.pyplot import *

class Bspline:
    def __init__(self, pts:'array(,2)', knots:'array(,1)'):
        self.pts = pts; self.knots = knots; self.n = len(knots) - 1
        self.d = [self.pts[0, 0], self.pts[-1, 0]]

    def getBasis(self, idx=0, deg = 1):
        '''Cox-de Boor formula'''
        x = self.knots.copy()

        def N(i, p) -> 'func':
            if p == 0:
                return lambda u: 1 if x[i] <= u and u <= x[i+1] else 0

            else:
                a1 = lambda u: u - x[i]; b1 = lambda _: x[i+p] - x[i]
                a2 = lambda u: x[i+p+1] - u; b2 = lambda _: x[i+p+1] - x[i+1]

                c1 = lambda u: 0 if abs(b1(u)) < 0.00001 else a1(u)/b1(u)
                c2 = lambda u: 0 if abs(b2(u)) < 0.00001 else a2(u)/b2(u)

                return lambda u: c1(u)*N(i, p-1)(u) + c2(u)*N(i+1, p-1)(u)

        return N(idx, deg)

    def renderBasis(self, deg) -> 'None':
        nsp = 1000
        N = iter(
                list(self.getBasis(i, deg)(u)
                     for u in linspace(self.d[0], self.d[1], nsp))
                 for i in range(self.n-deg))

        [plot(linspace(self.d[0], self.d[1], nsp), tmp) for tmp in N]

        xlim(-0.1, 1.1)
        ylim(-0.1, 1.1)
        grid()
        #show()



if __name__=='__main__':
    controlpts = array([[0, 0], [0.4, 0.2], [0.75, -0.2], [1,0]])
    knots = array([0,0, 0, 0.3, 0.5, 0.6, 1, 1, 1])
    B = Bspline(controlpts, knots)
    B.renderBasis(2)
    rcParams["figure.figsize"] = [16, 2]
    #savefig('basis3.pdf')
    show()
