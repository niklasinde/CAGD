from numpy import *
from matplotlib.pyplot import *

class Bspline:
    def __init__(self, pts:'array(,2)', knots:'array(,1)'):
        self.pts = pts; self.knots = knots

    def getBasis(self):
        '''Cox-de Boor formula'''
        x = self.knots.copy()

        def N(i, p):
            if p == 0:
                return lambda u: 1 if x[i] <= u and u <= x[i+1] else 0

            else:
                def c1(u):
                    if abs(x[i+p] - x[i]) <= 0.001: return 0
                    return (u - x[i])/(x[i+p] - x[i])
                def c2(u):
                    if abs(x[i+p+1] - x[i+1]) <= 0.001: return 0
                    return (x[i+p+1] - u)/(x[i+p+1] - x[i+1])

                return lambda u: c1(u)*N(i, p-1)(u) + c2(u)*N(i+1, p-1)(u)

        return N(0, len(self.knots)-2)

    def renderBasis(self):
        N = self.getBasis()
        sample = linspace(0, 1, 100)

        plot(sample, [N(_) for _ in sample])
        xlim(-0.1, 1.1)
        ylim(-0.1, 1.1)
        grid()
        show()



if __name__=='__main__':
    pass
