from numpy import *
from matplotlib.pyplot import *

class Bspline:
    def __init__(self, pts:'array(,2)', knots:'array(,1)'):
        self.pts = pts; self.knots = knots; self.n = len(knots) - 1

    def getBasis(self, xi:'float'):
        '''Cox-de Boor formula'''
        x = self.knots.copy()

        def N(i, p) -> 'func':
            if p == 0:
                return lambda u: 1 if x[i] <= xi and xi <= x[i+1] else 0

            else:
                a1 = lambda _: _ - x[i]; a2 = lambda _: x[i+p] - _
                b1 = lambda _: x[i+p-1] - x[i]; b2 = lambda _: x[i+p] - x[i+1]

                c1 = lambda u: a1(u)/b1(u) if abs(b1(u)) > 0.001 else 0
                c2 = lambda u: a2(u)/b2(u) if abs(b2(u)) > 0.001 else 0

                return lambda u: c1(u)*N(i, p-1)(u) + c2(u)*N(i+1, p-1)(u)

        return N(0, self.n-1)

    def renderBasis(self, SHOW=True) -> 'None':
        N = [self.getBasis(pt[0]) for pt in self.pts]
        sample = linspace(self.knots[0], self.knots[-1], 100)

        for b in N:
            plot(sample, [b(_) for _ in sample])

        #xlim(-0.1, 1.1)
        #ylim(-0.1, 1.1)
        grid()
        if SHOW == True: show()



if __name__=='__main__':
    controlpts = array([[-1, 0], [0.4, 0.2], [0.75, -0.2], [1,0]])
    knots = array([0, 0.25, 0.75, 1])
    B = Bspline(controlpts, knots)
    B.renderBasis(SHOW = True)
