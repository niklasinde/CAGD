#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:40:31 2017

@author: niklasinde
"""

from BernPoly import *
from BernRec import *





if __name__=='__main__':
    x = array([[0,0],
               [0.2, 0.5],
               [0.4, 0.2],
               [0.7, -0.2],
               [1, 0.1]])
    '''task1'''
    P = BernPoly(x)
#    P.renderBasisFuncs(0.2)
    '''task2'''
    pts3 = array([[0.05, 0.02], [0.1, 0.2],
                 [0.2, -0.1], [0.3, 0],
                 [0.4, 0.1], [0.7, 0.2]])
    #@timeit
    #def RecursiveMethod():
    #    proj1 = RecurBasis(pts)
    #    proj1.render()
    #@timeit
    #def BernsteinMethod():
    #    proj2 = BernPoly(pts)
    #    proj2.render()

    #RecursiveMethod()
    #BernsteinMethod()
    '''task3'''
     
    #B = BernPoly(pts)
    #B.render(domain=[-1,2], env=True)
    pts = array([[ -1.   ,  0.   ],[ 0.25 ,  1.   ],[ 0.375, -2.   ], [ 0.5 , -0.   ],[2,10]])
    A = BernPoly(pts)
    A.render()
    show()
    A.renderBasisFuncs(pt=True)
#    A.subdivision(t=0.4)

    #xlim(-0.2, 1.4)
    #ylim(-0.2, 0.4)
    grid()
    #savefig('SubDivision.pdf')
    show()