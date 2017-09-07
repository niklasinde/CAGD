from matplotlib.pyplot import *
from numpy import *

class casteljau():
    def __init__(self, pts, *kwargs):
        self.pts = pts

    def getB(self, i, k):
        if k == 0:
            return lambda t: self.pts[i-1, :]
        else:
            return lambda t: (1-t)*self.getB(i+1, k-1)(t) + t*self.getB(i, k-1)(t)

    def run(self, domain:'list [a, b]', nsp = 100) -> 'plot':
        b      = self.pts
        B      = self.getB(0, len(self.pts))
        sample = linspace(min(b[:, 0]), max(b[:, 0]), nsp)
        coords = array([list(B(_)) for _ in sample])

        plot(coords[:, 0], coords[:, 1], c='r')
        scatter(b[:, 0], b[:, 1])
        for _ in range(len(self.pts)-1):
            plot([b[_, 0], b[_+1, 0]], [b[_, 1], b[_+1, 1]], c='k')

        grid()

if __name__=='__main__':
    pts = array([[0.05, 0.05], [0.1, 0.2], [0.5, 0.3], [0.6, 0.2], [0.5, -0.1]])
    p = casteljau(pts)
    p.run([0, 0.6])
    show()





