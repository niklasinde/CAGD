from numpy import*
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

class Bezier:
    def __init__(self, points):
        self.b, self.n = points, len(points)-1
        self.domain = (min(points[:,0]), max(points[:, 0]))

    def __run(self):
        c1 = lambda t: (self.domain[-1] - t)/self.domain[-1] if self.domain[-1] != 0 else 0
        c2 = lambda t: (t - self.domain[0])/(self.domain[-1] - self.domain[0]) if (self.domain[-1] - self.domain[0]) != 0 else 0

        def B(i, p):
            if p == 0: return lambda _: self.b[i, :]
            return lambda t: c1(t)*B(i, p-1)(t) + c2(t)*B(i+1, p-1)(t)

        return B(0, self.n)

    def evaluate(self, nsp=100):
        basis  = self.__run()
        sample = linspace(self.domain[0], self.domain[-1], nsp)

        return array(list( basis(x) for x in sample ))

    def render(self):
        Y = self.evaluate()
        Axes3D.scatter(xs = Y[:, 0], ys = Y[:, 1], zs = Y[:,2])
        #for arg in self.b: Axes3D.scatter(*arg, c='k')
        #Axes3D.plot(self.b[:, 0], self.b[:, 1], self.b[:,2], alpha=0.2, c='k')


class surface:
    def __init__(self, points1, points2, points3, points4):
        self.b1, self.b2, self.b3, self.b4 = points1, points2, points3, points4



p1 = array([ [0,0,0], [1,2,1], [2,-1,1], [3,0,0] ])

C1 = Bezier(p1)
print(C1.evaluate())
C1.render()
plt.show()
