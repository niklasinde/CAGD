from numpy import*
from numpy.linalg import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


class BilinearPatch:
    def __init__(self, data, domain=(0,1)):
        self.data = data
        b = data
        print(b[0,0])
        print(b[0,1])
        print(b[1,0])
        print(b[1,1])

    def BilinearInterpolation(self):
        b = self.data

        def S(u, v):
            return (1-u)*(1-v)*b[0,0]\
                   + (1-v)*u*b[0,1]\
                   + u*(1-v)*b[1,0]\
                   + u*v*b[1,1]

        return S

    def evaluate(self):
        sample = linspace(0, 1, 10)
        S      = self.BilinearInterpolation()

        return (array(list( list( S(u, v) for u in sample) for v in sample )),
                array(list( list( S(v, u) for u in sample) for v in sample )) )


    def render(self):
        XY, YX = self.evaluate()

        for values in XY:
            print(values)

        for arg in self.data: ax.scatter(*arg, c='k')

if __name__=='__main__':
    p1 = array([ [-2,-2,1],
                 [2,2,1],
                 [1,0,0],
                 [0,1,0],
                 [1/2, 1/2, 1] ])

    domain = (0, 1)

    C = BilinearPatch(p1, domain)
    C.render()
    plt.show()
