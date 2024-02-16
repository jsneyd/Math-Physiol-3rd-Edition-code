
#   -------------------------------------------------------------------
# 
#    ODE integrator for the Noble Purkinje fiber model equations.
# 
#    For Chapter 12, Section 12.2.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# ODE RHS
def NobledeRHS(s, t):
    V, m, h, n = s   
    alphabeta = ab(V)
    currents = IV(s)
    INa = currents[0]
    IK = currents[1]
    Ian = currents[2]    
    Fv = -(INa + IK + Ian + Inp) / cm
    Fg = alphabeta[::2] * (1 - s[1:4]) - alphabeta[1::2] * s[1:4]    
    return [Fv, *Fg]

def ab(V):  
    # there are 5 parameters C for each alpha, beta
    alphabeta = (C[:,0]*np.exp(C[:,1]*(V-V0))+C[:,2]*(V-V0))/(1+C[:,3]*np.exp(C[:,4]*(V-V0)))
    return alphabeta


def IV(s):
    V, m, h, n = s   
    gna = gnabar * m ** 3 * h
    gk1 = 1.2 * np.exp(-(V + 90) / 50) + 0.015 * np.exp((V + 90) / 60)
    gk2 = 1.2 * n ** 4
    INa = gna * (V - Vna)
    Ik = (gk1 + gk2) * (V - Vk)
    Ian = gan * (V - Van)   
    return np.array([INa, Ik, Ian])


# Parameters
gnabar = 400.
Vna = 40.
gan = 0.14
Van = 40.
Vk = -100.
cm = 12.
# The matrix C contains parameters for the alpha's and beta's
C = np.array([[0,0,0.1,-1,-1/15],
              [0,0,-0.12,-1,0.2],
              [0.17,-0.05,0,0,0],
              [1,0,0,1,-0.1],
              [0,0,0.0001,-1,-0.1],
              [0.002,-0.0125,0,0,0]])
V0 = np.array([-48,-8,-90,-42,-50,-90])

tstep = 0.2
t_end = 3500

# Initial conditions
Inp = 0
V = -49
m = 0.1
h = 0.7
n = 0.3

# Time points
tspan = np.arange(0, t_end, tstep)

# Calculate the solution
s0 = [V, m, h, n]
S = odeint(NobledeRHS, s0, tspan)

# Plotting
T = tspan
V = S[:, 0]
m = S[:, 1]
h = S[:, 2]
n = S[:, 3]

plt.figure(1)
plt.subplot(2, 1, 1)
plt.plot(T, V, linewidth=2)
plt.xlabel('t (ms)', fontsize=20)
plt.ylabel('V (mV)', fontsize=20)

plt.subplot(2, 1, 2)
plt.plot(T, m, 'r', T, h, 'b', T, n, 'g', linewidth=2)
plt.legend(['m', 'h', 'n'], fontsize=16)
plt.xlabel('t (ms)', fontsize=20)

INa, IK, Icl = np.zeros_like(T), np.zeros_like(T), np.zeros_like(T)
for j, t in enumerate(T):
    s = S[j, :]
    currents = IV(s)
    INa[j] = currents[0]
    IK[j] = currents[1]
    Icl[j] = currents[2]

plt.figure(3)
plt.plot(T, INa, T, IK, T, Icl, linewidth=2)
plt.legend(['INa', 'IK', 'Il'], fontsize=16)
plt.xlabel('t (ms)', fontsize=20)


