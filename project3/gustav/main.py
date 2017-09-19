from numpy import *
from matplotlib.pyplot import *

from Basis import *

if __name__=='__main__':
    controlpts = array([[-1, 0], [0.4, 0.2], [0.75, -0.2], [1, 0]])
    knots = array([-1, 0, 1, 2])
    #knots = array([0, 0.25, 0.75, 1])
    B = Bspline(controlpts, knots)
    #B.getBasis(0)
    show()
    B.renderBasis(SHOW=True)

    #savefig('BasisTask1.pdf')

