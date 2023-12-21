#   -------------------------------------------------------------------
# 
#   This is an ode integrator for the HH equations.
# 
#    For Chapter 5, Fig. 5.6 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def deriv(t,y):
    v,n,m,h = y
    vd=(1.0/cm)*(-gk*(n**4)*(v-vk) - gna*(m**3)*h*(v-vna) - gl*(v-vl) + Iapp)
    md=alpham(v)*(1.0-m) - betam(v)*m
    hd=alphah(v)*(1.0-h) - betah(v)*h
    nd=alphan(v)*(1.0-n) - betan(v)*n
    return [vd,nd,md,hd]

def alpham(v):
    return 0.1*(25.0-v)/(np.exp((25.0-v)/10.0)-1.0)
def betam(v):
    return 4.0*np.exp(-(v)/18.0)
def alphan(v):
    return 0.01*(10.0-v)/(np.exp((10.0-v)/10.0)-1.0)
def betan(v):
    return 0.125*np.exp(-(v)/80.0)
def alphah(v):
    return 0.07*np.exp(-(v)/20.0)
def betah(v):
    return 1.0/(np.exp((30.0-v)/10.0)+1.0)

# initial conditions
v0=7.0
n0=0.31768
m0=0.052934
h0=0.59611

# the parameters
Iapp=0.0
gna=120.0
gk=36.0
gl=0.3
vna=115.0
vk=-12.0
vl=10.6
cm=1.0

# solve the odes
y0 = [v0,n0,m0,h0] 
tout = np.linspace(0,50,500)
soln = solve_ivp(deriv,[0, 50],y0,method='Radau',t_eval=tout,rtol=10e-6,atol=10e-9)

# plot the solution
plt.plot(soln.t,soln.y[0])
plt.xlabel(r'time (ms)')
plt.ylabel(r'$V-V_{\rm eq}$')
