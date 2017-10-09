from numpy import *
from matplotlib.pyplot import *
rcParams['figure.figsize'] = 8, 4


class RationalCasteljau:
    def __init__(self, pts):
        self.b = pts; self.n = len(pts)-1
        self.domain = (min(pts[:, 0]), max(pts[:, 0]))

    def run(self):
        points, weights = self.b, self.weights

        c1 = lambda t: (self.domain[-1] - t)/self.domain[-1]

        c2 = lambda t: (t - self.domain[0])/(self.domain[-1]
                                             - self.domain[0])


        def w(i, p):
            if p == 0: return lambda _: weights[i]

            return lambda t: c1(t)*w(i, p-1)(t) + c2(t)*w(i+1, p-1)(t)

        def N(i:'iter', p:'degree') -> 'func':
            if p == 0: return lambda _: points[i, :]

            return lambda t: (c1(t)*w(i, p-1)(t)*N(i, p-1)(t)
                              + c2(t)*w(i+1, p-1)(t)*N(i+1, p-1)(t))\
                             /w(i, p)(t)

        return N(0, self.n)


    def eval(self, w, color='k', alpha=1, env=False):
        self.weights = w

        basis = self.run()
        sample = linspace(self.domain[0], self.domain[-1], 100)

        y = array(list(basis(x) for x in sample))
        plot(y[:, 0], y[:, 1], c=color, alpha=alpha)
        if env==True:
            for arg in self.b: scatter(*arg, c=color)
            plot(self.b[:, 0], self.b[:, 1], alpha = 0.2, c='k')


if __name__=='__main__':
    #pts = array([[0, 0], [0.2, 0.5], [0.6, -0.2], [1, 0]])
    pts = array([ [0,0], [4,3], [3,1], [5,1] ])
    weights = array([1, 2, 3, 4])
    C = RationalCasteljau(pts)
    C.eval(weights, color='k', env=True)

    for J in range(1,10):
        J = 1/J
        weights = array([1, 2, 3, J])
        A = RationalCasteljau(pts)
        A.eval(weights, color='r', alpha=0.3)
    for J in range(1, 10):
        J = -1/J
        weights = array([1, J, 3, 4])
        A = RationalCasteljau(pts)
        #A.eval(weights, color='b', alpha=0.3)

    grid()
    #savefig('prob1.pdf')
    show()

