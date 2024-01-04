#   -------------------------------------------------------------------
# 
#    Code for the Pedersen model of secretory granule exocytosis.
# 
#    For Chapter 8, Section 8.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Define the right-hand side of the ODEs
def rhs(t,x):

    h = x[:ng]
    F = x[ng]
    I = x[ng + 1]

    M =  M1 * glucose(t)**m  / ( KM**m  +  glucose(t)**m ) +  M0  

# Warning! The order of the variables in trapz is different from Matlab. 
# This is a total pain in the neck.
    int_val = np.trapz( h[g < glucose(t)], g[g < glucose(t)])
    intinf = np.trapz(h,g)

    out = np.zeros_like(x)
    out[: ng] =  pplus * I *  phi  - pminus * h - fplus * h * (glucose(t) > g)
    out[ng] =  fplus * int_val -  k * F
    out[ng + 1] = M -  r * I -  pplus * I +  pminus * intinf

    return out

def glucose(t):
    if stim_choose == 1:
        G = [0, 50, 100, 150, 200]  # glucose levels
        times = [0, 2, 7, 12, 17]  # times of step application
        out = G[0] * (np.heaviside(t - times[0], 0) - np.heaviside(t - times[1], 0)) + \
              G[1] * (np.heaviside(t - times[1], 0) - np.heaviside(t - times[2], 0)) + \
              G[2] * (np.heaviside(t - times[2], 0) - np.heaviside(t - times[3], 0)) + \
              G[3] * (np.heaviside(t - times[3], 0) - np.heaviside(t - times[4], 0)) + \
              G[4] * np.heaviside(t - times[4], 0)

    if stim_choose == 2:
        G = [0, 300, 0, 300]  # glucose levels
        times = [0, 5, 60, 65]  # times of step application
        out = G[0] * (np.heaviside(t - times[0], 0) - np.heaviside(t - times[1], 0)) + \
              G[1] * (np.heaviside(t - times[1], 0) - np.heaviside(t - times[2], 0)) + \
              G[2] * (np.heaviside(t - times[2], 0) - np.heaviside(t - times[3], 0)) + \
              G[3] * np.heaviside(t - times[3], 0)

    return out


# Set initial conditions and parameters
Xmax   = 1.65
n   = 3.3
K   = 150
M1   = 1.25
m   = 12
KM   = 185
M0   = 0.02
fplus   = 6.2
k   = 1.09
pplus   = 0.021
pminus   = 0.11
r   = 0.0023
ng   = 600  # Number of values of g
g   = np.linspace(0, 600,  ng)
phi   =  n*g**(n - 1) *  K**n   / (( K**n   +  g**n)**2)
h0   =  Xmax   *  phi 

# Staircase experiment
init = np.concatenate([ [h0][0],[0,8] ])
stim_choose = 1  # 1 for staircase, 2 for repeated stimulus
tend = 22
tout = np.linspace(0,tend,1000)
sol = solve_ivp(rhs, [0,tend], init, method='BDF', t_eval=tout)

# Plot results
fig, ax1 = plt.subplots()
ax1.plot(sol.t,sol.y[ng],label='insulin')
plt.legend()
plt.xlabel('time(min)')
plt.ylabel('insulin secretion (\mug/min)')
ax2 = ax1.twinx()
ax2.plot(tout,glucose(tout),'red',label='glucose')
plt.ylabel('glucose (mg/100 ml)')
plt.legend()

# same stimulation twice
stim_choose = 2  # 1 for staircase, 2 for repeated stimulus
tend = 75
tout= np.linspace(0,tend,500);
sol = solve_ivp(rhs, [0,tend], init, method='LSODA', t_eval=tout,rtol=1e-10,atol=1e-10);

# Plot results
fig, ax1 = plt.subplots()
ax1.plot(sol.t,sol.y[ng])
plt.xlabel('time(min)')
plt.ylabel('insulin secretion (\mug/min)')
ax2 = ax1.twinx()
ax2.plot(tout,glucose(tout),'red')
plt.ylabel('glucose (mg/100 ml)')



