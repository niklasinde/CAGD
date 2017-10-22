from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


class deBoor:
    def __init__(self, points, knots, knots2, degree, degree2):
        self.b, self.u, self.v, self.p, self.q = points, knots, knots2,\
                                                 degree, degree2
        self.domain = [min(self.b[:, 0]), max(self.b[:, 0])]
        self.m = len(knots) - 1
        self.n = len(knots2) - 1

    def activeknots(self, x):
        for i in range(self.m):
            if self.u[i] <= x <= self.u[i+1]: return i
        else: raise ValueError(x, ' is not in the interval')

    def run(self, k:'position', x:'variable',
            t:'knots', c:'ctrlpts', p:'degree'):

        d = list(c[j+k-p] for j in range(0, p+1))

        def update(i, j):
            alpha = 1. if t[j+1+k-1] - t[j+k-p] == 0 else\
                    (x - t[j+k-p])/(t[j+1+k-i] - t[j+k-p])

            d[j] = (1 - alpha)*d[j-1] + alpha*d[j]

        list( list(update(i, j) for j in range(p, i-1, -1))
                                for i in range(1, p+1) )

        return d[p]

    def evaluate(self, x, controlpoints):
        if shape(array(x)) != ():
            return array(list((self.evaluate(_, controlpoints) for _ in x)))

        return self.run(self.activeknots(x), x,
                        self.u, controlpoints, self.p)

    def evaluate3D(self, x, y, controlpoints):
        if shape(array(y)) != ():
            innerpoints = self.evaluate(x, controlpoints)
            print(controlpoints)
            print(innerpoints)
            return array(list( (self.evaluate3D(x, _, innerpoints[i])
                               for i,_ in enumerate(y)) ))

        return self.run(self.activeknots(y), y,
                        self.v, controlpoints, self.q)

    def NURBS(self, weights):
        def step1():
            return array(list( weights[i]*append(x, 1)
                               for i,x in enumerate(self.b.copy()) ))

        def step2(i):
            X = linspace(self.u[i], self.u[i+1], 30)
            return self.evaluate(X, controlpoints = new_b)

        def step3(fx):
            return array(list( y[:-1]*y[-1] for y in fx ))

        new_b = step1()

        for i in range(self.m):
            Y = step3( step2(i) )
            ax.plot(Y[:, 0], Y[:, 1], Y[:, 2])


    def render(self,env = False,  alpha = 1):
        def PLOT(i):
            X = linspace(self.u[i], self.u[i+1], 30)
            Y = self.evaluate3D(X, X, controlpoints = self.b)
            print(Y)

            ax.plot(Y[:, 0], Y[:, 1], Y[:, 2], alpha =alpha)

        if env == True: list( scatter(*arg, c='k') for arg in self.b )
        list( PLOT(i) for i in range(self.m) )


if __name__=='__main__':
    points = array([ [0,0,0],
                     [0,0,0.5],
                     [0,0,1],
                     [0,0,2] ])

    knots = array([ 0, 0, 0, 0.5, 9, 9, 9])
    weights = array([ 1, 1, 1, 1])
    degree = 2
    A = deBoor(points, knots, knots, degree, degree)
    #A.NURBS(weights)
    A.render()

    grid()
    show()
