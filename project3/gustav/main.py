from numpy import *
from matplotlib.pyplot import *

from Basis import *

if __name__=='__main__':
    controlpts = array([[-1, 0], [0.4, 0.2], [0.75, -0.2], [1, 0]])
    knots = array([0, 0, 1, 1])
    B = Bspline(controlpts, knots)
    B.getBasis()

    #savefig('BasisTask1.pdf')

