#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def deriv(t,y):
    D,N,v = y
    fD = kD0 + kkD*v**2/(v**2+Kv1**2)
    fN = kN0 + kkN*(v**2/(v**2+Kv2**2))*(KD**2/(KD**2 + D**2))
    I = I0 + pv*v*N**2/(N**2 + KN**2)
    fv = v0*(1-np.exp(-I))/(1+np.exp(-I))
    Ddiff = fD - dD*D
    Ndiff = fN - dN*N
    vdiff = fv - dv*v
    return [Ddiff,Ndiff,vdiff]


I0 = 0.08
dD = 0.25
dN = 1
dv = 10
kkD = 4.5
kkN = 320
kD0 = 0.175
kN0 = 0
pv = 1
v0 = 30000
KD = 0.3
KN = 32
Kv1 = 1200
Kv2 = 1200

y0 = [0.7,0,0] 
tout = np.linspace(0,140,10000)
# solve the odes
soln = solve_ivp(deriv,[0, 140],y0,method='Radau',t_eval=tout,rtol=10e-6,atol=10e-9)
plt.plot(soln.t,soln.y[2])

np.savetxt('test.dat', np.transpose([soln.t,soln.y[0],soln.y[2]]) )
