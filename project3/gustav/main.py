from numpy import *

from Basis import *

if __name__=='__main__':
    controlpts = array([[0, 0], [0.5, 0.2], [0.75, -0.2], [1, 0]])
    knots = array([0, 0,0, 1, 1, 1])
    B = Bspline(controlpts, knots)
    B.getBasis()
    B.renderBasis()

