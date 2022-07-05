#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def deriv(t,y):
    v,g,z = y
    h = np.maximum(v,0)
    p = p1/(1+np.exp(-p2*(h-p3)))
    vd = s0*(-v*(v-c)*(v-1) - k1*g + k2*a(t))
    gd = k3*a(t) + k4*v - rg*g
    zd = p - rz*z
    return [vd,gd,zd]

def a(t):
    out=0
    if (np.remainder(t,stimperiod)<stimwidth):
        out=4
    return out

    
s0=50
c=0.2
k1=1
k2=0.02
k3=0.02
k4=1
rg=2.5
p1=100
p2=100
p3=0.3
rz=1

stimperiod = 1
stimwidth = 0.9*stimperiod

y0 = [0.1,0,0] 
tout = np.linspace(0,15,100000)
# solve the odes
soln = solve_ivp(deriv,[0, 15],y0,method='Radau',t_eval=tout,rtol=10e-6,atol=10e-9)
plt.plot(soln.t,soln.y[2])
