#     -------------------------------------------------------------------
# 
#    Model of pulsatile GnRH secretion by Brown et al.
# 
#      For Chapter 16, Section 16.2.2 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

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
        out=stimheight
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

stimperiod = 0.6
stimwidth = 0.9*stimperiod
stimheight = 4

y0 = [0.1,0,0] 
tout = np.linspace(0,15,1000)
# solve the odes
soln = solve_ivp(deriv,[0, 15],y0,method='BDF',t_eval=tout)
plt.plot(soln.t,soln.y[2])
