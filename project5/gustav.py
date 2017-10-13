%matplotlib inline
from numpy import *
from matplotlib.pyplot import *
from numpy.linalg import solve

class Bspline:
    def __init__(self, points, knots, degree):
        self.b, self.u, self.p, self.n = points, knots, degree, len(points)
        
    def run(self) -> 'func':
        u = self.u.copy()
        
        def div(lhs, rhs):
            if rhs == 0: return 0
            
            return lhs / rhs
        
        def recursion(i, k, t):
            if k == 0:
                if u[i] == u[i+1]: return 0
                if not (u[i] <= t <= u[i+1]): return 0
                
                return 1
            
            return div( (t - u[i]), (u[i+k] - u[i]) )*recursion(i, k-1, t)\
                   + div( (u[i+k+1] - t), (u[i+k+1] - u[i+1]) )*recursion(i+1, k-1, t)
        
        def N(i, k) -> 'func':
            return lambda t: recursion(i, k, t)
        
        return N
    
    def coefficients(self, interpolate):
        if interpolate == False: return self.b
        
        u = self.u.copy()
        N = self.run()
        M = zeros( (self.n, self.n) )
        
        def addBasis(row, col, value):
            M[row, col] = value
            
        def parameter(col):
            return u[0] + col*(u[-1] - u[0])/(self.n-1)
        
        [[addBasis(row, col, N(col, self.p)(parameter(row))) for row in range(self.n)] for col in range(self.n)]
        
        X, Y = solve(M, self.b[:, 0]), solve(M, self.b[:, 1])
        return array(list( [X[i], Y[i]] for i in range(self.n) ))

    def evaluate(self, interpolate=False, t=False, nsp=100):
        b = self.b.copy()
        N = self.run()
        
        basisvector = array(list( N(i, self.p) for i in range(len(self.u) - self.p - 1) ))
        coeffs = self.coefficients(interpolate)
        #for args in coeffs: scatter(*args, c='r')
        
        getValue = lambda t: sum( coeffs[i]*basis(t)
                                  for i,basis in enumerate(basisvector) )
        
        if t != False: return getValue(t)
        return array(list( getValue(t) for t in linspace(self.u[0], self.u[-1], nsp) ))
        
        

#for args in points: scatter(*args, c='k')
#plot(Y[:, 0], Y[:, 1], c='k')

import time as time

t0 = time.time()
knots = np.array([0, 0, 0, 3, 6, 6, 6])
points = np.array([ [0, 0], [6, 10], [7, 10.2], [9, 8] ])
C = Bspline(points, knots, 2)
Y = C.evaluate(interpolate=True)
plot(Y[:, 0], Y[:, 1], c='k')
show()
print(time.time()-t0,"time")

