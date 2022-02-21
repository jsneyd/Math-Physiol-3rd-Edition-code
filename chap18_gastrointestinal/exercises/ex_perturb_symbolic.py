#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 14:55:50 2022

@author: james
"""

from sympy import *

m0,a,y,f1,f2,beta = symbols('m0 a y f1 f2 beta')

f1 = -m0*y*y + beta + m0*a*(2-a)
f2 = -2*m0*a*y + beta + 2*m0*a

print(simplify(integrate(f1,(y,0,a)) + integrate(f2,(y,a,1))))
