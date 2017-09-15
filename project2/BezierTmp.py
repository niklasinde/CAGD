from numpy import *
from scipy.special import binom
from matplotlib.pyplot import *
from functools import reduce
import time

def timeit(method):
    def timed(*args, **kw):
        ts     = time.time()
        result = method(*args, **kw)
        te     = time.time()

        print('%r (%r, %r) %2.2f sec' % (method.__name__, args, kw, te-ts))
        return result
    return timed


class BernBasis:
    """
    Generates one of the Bernoubil basis. Gets called by the BernPoly in the function
    _generateBasis.
    """
    def __init__(self, degree:'int', iteration:'int'):
        self.n, self.i = degree, iteration
        assert self.n >= self.i, 'passed max iter.'

    def __call__(self, arg:'float') -> 'float':
        var = lambda x: x**self.i * (1 - x)**(self.n - self.i)
        return binom(self.n, self.i)*var(arg)
class RecurBasis:
    '''Used for task 3, to mimic the
       computation of homework 1.'''
    def __init__(self, pts):
        self.pts, self.n = pts, len(pts)-1
        self.Bezier = self._generateCurve()

    def _generateCurve(self) -> 'func':
        def b(i, k):
            if k == 0: return lambda _: self.pts[i, :]
            else: return lambda t: (1-t)*b(i, k-1)(t) + t*b(i+1, k-1)(t)

        return b(0, self.n)

    def _render(self, nsp=100, colour='r', env=None) -> 'None':
        values = array(list(self.Bezier(x) for x in linspace(0, 1, nsp)))

        plot(values[:, 0], values[:, 1], c=colour)
        if env:
            scatter(self.pts[:, 0], x[:, 1], c='k', alpha=0.5)
            plot(list([self.pts[i, 0], self.pts[i+1, 0]] for i in range(self.n)),
                 list([self.pts[i, 1], self.pts[i+1, 1]] for i in range(self.n)),
                 c='k', alpha=0.2)


class BernPoly:
    """
    Calls the "BernBasis class
    """
    def __init__(self, pts:'vector'):
        self.pts, self.n = pts, len(pts) - 1
        self.basis = self._generateBasis()
        self.px    = self._generatePolynomial()
        self.domain = [pts[0][0],pts[-1][0]]

    def __call__(self, x) -> 'matrix':
#        if not (0 <= x <= 1).any(): raise ValueError("x is not in the interval")
        return self.px(x)

    def _generatePolynomial(self) -> 'func':
        return lambda x: sum(self.pts[i]*base(x)
                             for i, base in enumerate(self.basis))

    def _render(self, nsp=100, color='r', env=True) -> 'None':
        domain = self.domain
        print(domain[0])
        sample = linspace(domain[0], domain[1], nsp)
        u = lambda a: a
        values = array(list(self.px(u(x))
                       for x in linspace(0, 1, nsp)))

        plot(sample, values[:, 1], c=color)

        if env:
            scatter(self.pts[:, 0], self.pts[:, 1], c='k', alpha=0.5)
            plot(list([self.pts[i, 0], self.pts[i+1, 0]] for i in range(self.n)),
                 list([self.pts[i, 1], self.pts[i+1, 1]] for i in range(self.n)),
                 c='k', alpha=0.2)

    def _PlotBasisFuncs(self, pt=False) -> 'None':
        tmp = self.basis.copy()
        sample = linspace(0, 1, 100)
        for i, p in enumerate(tmp):
            plot(sample, vectorize(p)(sample))

        xlim(-0.2, 1.2)
        ylim(-0.2, 1.2)

        if pt: print(sum(base(pt) for base in tmp))

    def _generateBasis(self) -> 'matrix':
        return array(list(BernBasis(self.n, idx) for idx in range(self.n + 1)))
   
    
    def addpoint(self, p:"vetor:point(s) to add",pos:"int possision of the point"):
        if np.shape(p)[1]==2: # Check correct shape
            raise ValueError("Input array is the wrong shape. Has to be np.array([[x,y],[x_1,y_1]....])")
        else:
            self.pts = np.insert(selt.pts,pos,inp,axis=0) #update basis and polynomial.
            self.basis = self._generateBasis()
            self.px    = self._generatePolynomial()
        return(self.pts)
    def movepoint(self, p:"vetor:point(s) to add",pos:"int possision of the point"):
        if np.shape(p)[1]==2: # Check correct shape
            raise ValueError("Input array is the wrong shape. Has to be np.array([[x,y],[x_1,y_1]....])")
        else:
            self.pts[pos] = p #update basis and polynomial.
            self.basis = self._generateBasis()
            self.px    = self._generatePolynomial()
        return(self.pts)
            
class Subdiv:
    def __init__(self, Bernpoly, t0):
        """
        Maybe we should make sub div as a subclass to Bernpoly?
        Inputs: Points info
        
        This functions create two new bernpoly with a start and end point in t0
        We have to implement this algoritm for it to work.
        https://math.stackexchange.com/questions/1408478/subdividing-a-b%C3%A9zier-curve-into-n-curves
        """
        self.t0 = t0
        self.bernorginal, self.pts = Bernpoly, Bernpoly.pts
        
        self.ht = self.hot_interval(self.pts,self.t0)
        self.newpoints()
#        print(self.pts1,"\n",self.pts2)
        self.updatebern()

    def __call__(self):
        return self.Bern1,self.Bern2
    
    def hot_interval(self,pts,val):
        pts = [x[0] for x in pts]
        lo, hi = 0, len(pts) - 1
        while lo <= hi:

            mid = (lo + hi) // 2
#            print(mid,pts[mid],pts,val)
            if pts[mid] == val == pts[mid+1]:
                return mid
            if pts[mid] < val < pts[mid+1]:
                return mid
            elif pts[mid] < val: 
                lo = mid + 1
            elif val < pts[mid]:
                hi = mid - 1
            else:
                print("Something is wrong")
                
    def newpoints(self):
        pts1,pts2 = self.pts[0:self.ht+1,:],self.pts[self.ht+1:,:] #Slice the original points in to two groups
        print(len(pts1))
        self.pts1, self.pts2 = self.fixpoints(pts1,np.shape(pts1)[0]),self.fixpoints(pts2,0)

    
    
    
    
    def fixpoints(self,points, pos):
        print(self.t0,self.bernorginal(self.t0))
        inp = np.array([self.bernorginal(self.t0)])
        pts = np.insert(points,pos,inp,axis=0) #update basis and polynomial.
        return(pts)
    def updatebern(self):
        self.Bern1 = BernPoly(self.pts1)
        self.Bern2 = BernPoly(self.pts2)
        
    def _render(self,env=True) -> 'None':
        self.Bern1._render(env=True)
        self.Bern2._render(env=True,color="g")
        show()
            


'''task1'''
x = array([[0,0], [0.2, 0.5], [0.4, 0.2], [0.4, 0.2], [0.7, -0.2], [1, 0.1]])


def task1():
    P = BernPoly(x)
    P._PlotBasisFuncs(0.2)
#task1()
'''task2'''
pts = array([[0, 0], [0.25, 0.5], [0.75, 0.5], [1, 0]])
def task2():
    @timeit
    def RecursiveMethod():
        proj1 = RecurBasis(pts)
        proj1._render()
    @timeit
    def BernsteinMethod():
        proj2 = BernPoly(pts)
        proj2._render()
    
    RecursiveMethod()
    BernsteinMethod()
#task2()
'''task3'''
def task3(save=False):
    A = BernPoly(pts)
    A._render(domain=[0,1], env=True)
    
    xlim(-0.2, 1.4)
    ylim(-0.2, 0.4)
    grid()
    if save:
        savefig('basisfuncs.pdf')
    show()
#task3()
    
'''Task4'''
#def task4():
Bern1 = BernPoly(pts)
Bern1._render()
split = 0.5
a = Subdiv(Bern1,split)
a._render(env=True)
bern1,bern2 = a()
#print(bern1.pts,"\n",bern2.pts)
#show()

#Bern2._render()
    
#task4()
    
