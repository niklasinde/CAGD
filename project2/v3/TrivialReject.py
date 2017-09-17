from numpy import *
from matplotlib.pyplot import *

class rectangle:
    def __init__(self, x, y):
        if abs(x[0] - x[1]) < 0.01: raise Exception
        self.x = x; self.y = y

    def checkBoundary(self, A, B):
        for a in A:
            return True if a in set(B) else False
        else: False

    def checkInside(self, A, B):
        return not (A[0:2] == B[0:2] or A[0:2] == B[2:4])

    def control(self, other:'obj'):
        X = list(self.x).copy(); [X.append(_) for _ in other.x]
        sX = X.copy(); sX.sort()

        Y = list(self.y).copy(); [Y.append(_) for _ in other.y]
        sY = Y.copy(); sY.sort()

        d = 0.2
        plot(X[0:2], [Y[0], Y[0]], c='k', alpha=d)
        plot(X[0:2], [Y[1], Y[1]], c='k', alpha=d)
        plot([X[0], X[0]], Y[0:2], c='k', alpha=d)
        plot([X[1], X[1]], Y[0:2], c='k', alpha=d)

        if self.checkBoundary(self.x, other.x) or self.checkInside(sX, X):
            if self.checkBoundary(self.y, other.y) or self.checkInside(Y, sY):
                return True
            else: return False
        else: return False


def drawCurves(X):
    plot(X[:, 0], X[:, 1])

if __name__=='__main__':
    pts = array([[0,0], [1,1]])
    ptsline = array([[-0.5, 0.5], [0.1, -0.1]])
    sq1 = rectangle(pts[:, 0], pts[:, 1])
    t = sq1.control(rectangle(ptsline[:, 0], ptsline[:, 1]))
    print(t)
    xlim(-1.5, 3.5); ylim(-0.5, 4.5)
    drawCurves(pts)
    drawCurves(ptsline)
    show()
