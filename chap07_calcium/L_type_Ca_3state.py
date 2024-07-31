#   -------------------------------------------------------------------
# 
#    3-state model of a L-type calcium channel
# 
#    For Chapter 7, Section 7.4.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def Ltrhs(t, x):
    c = 2 if t > 3 else 0
    gamma = 0.44 * c
    gammabar = gamma * (beta + 2 * alpha) / (alpha + beta)
    
    op = x[0]
    n3 = x[1]
    c3 = 1 - n3 - op
    
    opp = -g * op + fbar * n3
    n3p = -(fbar + gammabar) * n3 + omegabar * c3 + g * op
    
    return [opp, n3p]


f = 0.85
omega = 0.02
g = 2

Vmlist = [-75, -55, -35, -15, 5, 25, 45]
keep = []

plt.figure()
for Vm in Vmlist:
    alpha = 2 * np.exp(0.012 * (Vm - 35))
    beta = 0.0882 * np.exp(-0.05 * (Vm - 35))
    c = 0
    gamma = 0.44 * c
    
    fbar = alpha * f / (alpha + beta)
    omegabar = omega * (beta / 2 + alpha) / (2 * alpha + beta / 2)
    
    dt = 0.01
    tend = 10
    init = [0.0102, 0.9898]
    tspan = np.arange(0, tend + dt, dt)
    
    sol = solve_ivp(lambda t, x: Ltrhs(t, x), [0, tend], init, t_eval=tspan, method='LSODA')
    op = sol.y[0, :]
    keep.append(op)
    
    plt.plot(sol.t, op, label=f'V={Vm}')


plt.xlabel('t (ms)')
plt.ylabel('Open probability,o')
plt.legend()


