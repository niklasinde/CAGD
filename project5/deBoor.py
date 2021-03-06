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

        def step1():
            return  array(list( weights[i]*append(x, 1) for i,x in enumerate(self.b.copy()) ))

        def step2(i):
            X = linspace(self.u[i], self.u[i+1], 30)
            return self.eval(X, controlpoints = new_b)

        def step3(fx):
            return array(list( y[:2]*y[-1] for y in fx ))

        new_b = step1()

        for i in range(self.m):
            Y = step3( step2(i) )

            plot(Y[:, 0], Y[:, 1])


    def render(self,env = False,  alpha = 1):
        if env == True:
            [scatter(*arg, c='k') for arg in self.b]

        for i in range(self.m):
            X = linspace(self.u[i], self.u[i+1], 30)
            Y = self.eval(X, controlpoints = self.b)

            plot(Y[:, 0], Y[:, 1], alpha =alpha)

class parameters:
    def __init__(self, points):
        self.b = points[:, 0]
        self.n = len(points)

    def uniform(self):
        a, b = min(self.b), max(self.b)
        print(a, b, self.n)
        tmp = array([a])
        for k in range(1, self.n):
            midpoints = lambda j: (a + j*(b - a)/(self.n - 1))
            _tmp = midpoints(k)
            print(_tmp)
            tmp = append(tmp, _tmp)
        tmp = append(tmp, b)
        return tmp

if __name__=='__main__':
    points = array([ [0,0],
                     [6,10],
                     [7, 10.2],
                     [9, 8] ])
    knots = array([ 0, 0, 0, 6, 9, 9, 9])
    degree = 2

    A = parameters(points)
    B = A.uniform()
    knots = array([0, 0])
    knots = append(knots, B)
    knots = append(knots, 9)
    print(knots)
    print(B)

    Curve = deBoor(points, knots, degree)
    Curve.render(env = True, alpha=0.5)

    grid()
    #savefig('prob4.pdf')
    show()
