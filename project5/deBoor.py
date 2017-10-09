from numpy import *
from matplotlib.pyplot import *


class deBoor:
    def __init__(self, points, knots, degree):
        self.b, self.u, self.p = points, knots, degree
        self.domain = [min(self.b[:, 0]), max(self.b[:, 0])]
        self.m = len(knots) - 1

    def activeknots(self, x):
        for i in range(self.m):
            if self.u[i] <= x <= self.u[i+1]: return i
        else: raise ValueError(x, ' is not in the interval')

    def run(self, k, x, t, c, p):
        d = list(c[j+k-p] for j in range(0, p+1))

        for i in range(1, p+1):
            for j in range(p, i - 1, -1):
                if t[j+1+k-i] - t[j+k-p] == 0:
                    alpha = 1
                else:
                    alpha = (x - t[j+k-p]) / (t[j+1+k-i] - t[j+k-p])

                d[j] = (1 - alpha)*d[j-1] + alpha*d[j]

        return d[p]

    def eval(self, x, controlpoints):
        if shape(array(x)) != (): return array(list((self.eval(_, controlpoints) for _ in x)))

        return self.run(self.activeknots(x), x, self.u, controlpoints, self.p)

    def NURBS(self, weights):
        # (1) Multiply control points with weights.
        # (2) Run de Boor.
        # (3) Convert back by dividing the first components by the last one.
        # (-) The last component of each new control point is its weight.

        new_b = array(list( weights[i]*append(x, 1) for i,x in enumerate(self.b.copy()) ))

        for i in range(self.m):
            X = linspace(self.u[i], self.u[i+1], 30)
            Y = self.eval(X, controlpoints = new_b)

            Y = array(list( y[:2]*y[-1] for y in Y ))

            plot(Y[:, 0], Y[:, 1])


    def render(self):
        for i, arg in enumerate(self.b):
            if i in [1, 3, 5, 7]: scatter(*arg, c='r')
            else: scatter(*arg, c='k')


        for i in range(self.m):
            X = linspace(self.u[i], self.u[i+1], 30)
            Y = self.eval(X, controlpoints = self.b)

            plot(Y[:, 0], Y[:, 1], alpha =0.2)

if __name__=='__main__':
    points = array([ [0.5, 0],
                     [0, 0],
                     [0, 0.5],
                     [0, 1],
                     [0.5, 1],
                     [1, 1],
                     [1, 0.5],
                     [1, 0],
                     [0.5, 0] ])

    knots = array([ 0, 0, 0, 0.15, 0.35, 0.50, 0.70, 0.85, 0.95, 1, 1, 1])
    weights = array([ 1,
                      0.9,
                      0.65,
                      0.88,
                      1.03,
                      0.92,
                      1.01,
                      0.95,
                      1 ])
    degree = 2

    print('nr of knots: ', len(knots), 'nr of points: ', len(points[:, 0]), 'degree: ', degree,
            'm: ', len(points[:, 0])+degree, '|||||| nr of weights :', len(weights))

    Curve = deBoor(points, knots, degree)
    Curve.render()
    Curve.NURBS(weights)
    grid()
    #savefig('prob3_1.pdf')
    show()



    '''
    points2 = array([ [0.7, -0.4], [1, -0.4], [2.5, -1.2],
                     [3.2, -0.5], [-0.2, -0.5], [0.5, -1.2],
                     [2, -0.4], [2.3, -0.4] ])

    knots2 = array([ 1, 1, 1, 1, 6/5, 7/5, 8/5, 9/5, 2, 2, 2, 2])
    '''
