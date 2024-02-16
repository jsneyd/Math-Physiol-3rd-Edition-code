
#   -------------------------------------------------------------------
# 
#    ODE integrator for the Beeler-Reuter equations.
# 
#    For Chapter 12, Section 12.2.3 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def BRode(t, s):    
    V, m, h, j, d, f, x, Cai = s
    gate = [m,h,j,d,f,x]
    
    # get all the alpha's and beta's
    alpbet = ab(V)
    A = alpbet[0:11:2]
    B = alpbet[1:12:2]
    
    Inf = A / (A + B)
    Tau = 1 / (A + B) 
    C_m = 1
    BCL = 1000
    I = 10 if 10 < (t % BCL) < 15 else 0
    
    currents = IV(s)
    INa, IK, Ix, Is = currents

    V_prime = np.array([-(INa + IK + Ix + Is - I) / C_m])
    gate_prime = A * (np.ones(6) - gate) - B * gate
    Cai_prime = np.array([0.07 * (10**(-7) - Cai) - 10**(-7) * Is])
    
    s_prime = np.concatenate((V_prime, gate_prime, Cai_prime))
    return s_prime

def IV(s):
    V, m, h, j, d, f, x, Cai = s
    
    g_Na = gNa0 * m**3 * h * j + 0.003
    INa = g_Na * (V - ENa)
    Ek1 = EK
    Ix = 0.8 * x * (np.exp(0.04 * (V + 77)) - 1) / np.exp(0.04 * (V + 35))
    IK1 = 1.4 * (np.exp(0.04 * (V + 85)) - 1) / (np.exp(0.08 * (V + 53)) + np.exp(0.04 * (V + 53))) + 0.07 * (V + 23) / (1 - np.exp(-0.04 * (V + 23)))
    IK = IK1
    Esi = -82.3 - 13.0287 * np.log(Cai)
    Is = 0.09 * d * f * (V - Esi)
    
    currents = [INa, IK, Ix, Is]
    return currents

def ab(V):
    return (C[:, 0] * np.exp(C[:, 1] * (V - V0)) + C[:, 2] * (V - V0)) / (1 + C[:, 3] * np.exp(C[:, 4] * (V - V0)))

# Constants

RTbyF = 25.8
F = 96485.
CNai = 50
CNae = 457
CKi = 397
CKe = 20
ENa = RTbyF * np.log(CNae / CNai)
EK = RTbyF * np.log(CKe / CKi)
ENa = 50
gNa0 = 4
gK0 = 0.8
gK10 = 0.35
V0=[-47,-72,-77,-22.5,-78,-32,5,-44,-28,-30,-50,-20]

C = np.array([
    [0, 0, 1, -1, -0.1],
    [40, -0.056, 0, 0, 0],
    [0.126, -0.25, 0, 0, 0],
    [1.7, 0, 0, 1, -0.082],
    [0.055, -0.25, 0, 1, -0.2],
    [0.3, 0, 0, 1, -0.1],
    [0.095, -0.01, 0, 1, -0.072],
    [0.07, -0.017, 0, 1, 0.05],
    [0.012, -0.008, 0, 1, 0.15],
    [0.0065, -0.02, 0, 1, -0.2],
    [0.0005, 0.0833, 0, 1, 0.057],
    [0.0013, -0.06, 0, 1, -0.04]
])

# Steady state initial condition:
s0=[ -84.5738,
	0.011,
	0.9877 ,
	0.9748,
	0.003,
	1,
	0.0056,
	0.0000001782]

# Time parameters
t_end = 2000
tspan = np.linspace(0, t_end,10000)

# Solve the ODE system
# You need to be really careful with Python ODE solvers, and set the max_step small.
# Otherwise it will take huge steps and miss critical features.
sol = solve_ivp(BRode,[0,t_end],s0,method ='Radau',t_eval=tspan,max_step=1)

# Plot the results
fig, axs = plt.subplots(2, 1, figsize=(10, 8))
axs[0].plot(sol.t, sol.y[0], linewidth=2)
axs[0].set_xlabel('t (ms)', fontsize=20)
axs[0].set_ylabel('V (mV)', fontsize=20)
axs[1].plot(sol.t, sol.y[1],sol.t, sol.y[2],sol.t, sol.y[3],sol.t, sol.y[4],
            sol.t, sol.y[5],sol.t, sol.y[6],sol.t, sol.y[7],linewidth=2)
axs[1].legend(['m', 'h', 'j', 'd', 'f', 'x', 'Cai'], fontsize=16)
axs[1].set_xlabel('t (ms)', fontsize=20)
plt.show()

# Compute currents
currents = IV(sol.y)
INa, IK, Ix, Is = currents

# Plot the currents
fig, axs = plt.subplots(2,1)
axs[0].plot(sol.t, IK, 'r', sol.t, Ix, 'b', sol.t, Is, 'm', linewidth=2)
axs[0].legend(['IK', 'Ix', 'Is'], fontsize=16)
axs[0].set_xlabel('t (ms)', fontsize=20)
axs[0].set_ylabel('Current (nA)', fontsize=20)
axs[1].plot(sol.t, INa, 'k', linewidth=2)
axs[1].legend(['INa'], fontsize=16)
axs[1].set_xlabel('t (ms)', fontsize=20)
axs[1].set_ylabel('Current (nA)', fontsize=20)
plt.show()
