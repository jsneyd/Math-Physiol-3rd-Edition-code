#   ------------------------------------------------------------------- 
#    This is a simple ode integrator for the reduced HH equations (also called
#    the fast-slow HH equations). 
#
#    For Chapter 5, Fig. 5.10 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def deriv(t,y):
    v,n = y
    minf = alpham(v)/(alpham(v)+betam(v))
    vd=(1.0/cm)*(-gk*(n**4)*(v-vk) - gna*(minf**3)*(0.8-n)*(v-vna) - gl*(v-vl) + Iapp)
    nd=alphan(v)*(1.0-n) - betan(v)*n
    return [vd,nd]

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

def nullcline(n,v):
    minf = alpham(v)/(alpham(v)+betam(v))
    return (1.0/cm)*(-gk*(n**4)*(v-vk) - gna*(minf**3)*(0.8-n)*(v-vna) - gl*(v-vl) + Iapp)

v0=15
n0=0.31768
m0=alpham(v0)/(alpham(v0)+betam(v0))
h0=0.8-n0

Iapp=0.0
gna=120.0
gk=36.0
gl=0.3
vna=115.0
vk=-12.0
vl=10.6
cm=1.0

# Solve the odes and plot the solution
y0 = [v0,n0] 
tout = np.linspace(0,30,5000)
soln = solve_ivp(deriv,[0, 30],y0,method='Radau',t_eval=tout,rtol=10e-6,atol=10e-9)
plt.plot(soln.t,soln.y[0])

# Calculate the nullclines
vv = np.linspace(-11,112,100)
m_nullcline = alphan(vv)/(alphan(vv)+betan(vv))
n_nullcline = np.zeros(100)
mm = np.zeros(100)
for i in range(100):
    initial_guess = 0.5
    n_nullcline[i] = fsolve(nullcline, initial_guess,args=(vv[i],))

# plot the nullclines and add the solution in the phase plane
plt.figure()
plt.plot(vv,n_nullcline,'--g',label='dv/dt=0')
plt.plot(vv,m_nullcline,'--b',label='dn/dt=0')
plt.plot(soln.y[0],soln.y[1],'r',label='solution')
plt.xlabel('v (mV)')
plt.ylabel('n')
plt.legend()