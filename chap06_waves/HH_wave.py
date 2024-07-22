
#   -------------------------------------------------------------------
# 
#    Code to simulate HH eqns via Method of lines
# 
#    For Chapter 6, Figure 6.11 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def deRHS(t, s):
    V = s[:N]
    m = s[N:2 * N]
    h = s[2 * N:3 * N]
    n = s[3 * N:4 * N]

    AM = .1 * (25. - V) / (np.exp(.1 * (25. - V)) - 1.)
    BM = 4. * np.exp(-V / 18.0)
    AH = 0.07 * np.exp(-V / 20.)
    BH = 1.0 / (np.exp(0.1 * (30. - V)) + 1.)
    AN = 0.01 * (10.001 - V) / (np.exp(.1 * (10.001 - V)) - 1.)
    BN = 0.125 * np.exp(-V / 80.)

    Fm = (AM * (1. - m) - BM * m)
    Fh = (AH * (1. - h) - BH * h)
    Fn = (AN * (1. - n) - BN * n)

    INa, IK, ICl = IV(V, m, h, n)
    Ion = -INa - IK - ICl

    Fv = dg * (-sc * V + np.concatenate([[0], V[:-1]]) + np.concatenate([V[1:], [0]])) / taum + Ion / Cm
    return np.concatenate([Fv, Fm, Fh, Fn])

def IV(V, m, h, n):
    gna = gnabar * m ** 3 * h
    gk = gkbar * n ** 4
    ICl = gl * (V - Vl)
    INa = gna * (V - Vna)
    IK = gk * (V - Vk)
    return INa, IK, ICl


# Define global parameters
Veq = 0
gnabar = 120.  # mSiemens /cm^2
gkbar = 36.
gl = 0.3

Cm = 1  # mu F/cm^2
Tfact = 1  # temperature correction factor
taum = 1  # time constant in ms
lambda_m = 6.5  # space constant in mm.

Vna = Tfact * (115 + Veq) - Veq
Vk = Tfact * (-12 + Veq) - Veq
Vl = Tfact * (10.5988 + Veq) - Veq

L = 80
N = 101  # number of grid points
h = L / (N - 1)
dg = 1 / (h ** 2)  # coupling coefficient
X = np.linspace(0, L, N)

# Initial conditions
V0 = 20.1 * (1 - np.tanh(X / 2))
m0 = 0.0529
h0 = 0.5961
n0 = 0.3177
u0 = np.concatenate([V0, m0 * np.ones(N), h0 * np.ones(N), n0 * np.ones(N)])

sc = np.concatenate([[1], 2 * np.ones(N - 2), [1]])  # diagonal entries for the diffusion matrix
tstep = 0.05
t_end = 60
tspan = np.arange(0, t_end + tstep, tstep)

sol = solve_ivp(deRHS, [0, t_end], u0, t_eval=tspan, method='LSODA', max_step=1)

# Plotting
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
X_mesh, T_mesh = np.meshgrid(X, sol.t)
ax.plot_surface(X_mesh, T_mesh, np.transpose(sol.y[:N,:]), cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('T')
ax.set_zlabel('V')

# Now find the speed
thresh = 40
Tc = np.array([np.argmax(sol.y[j,:] >= thresh) * tstep for j in range(N)])
p = np.polyfit(Tc, X, 1)
speedest = p[0]
spest = p[1] + p[0] * Tc

plt.figure(2)
plt.plot(Tc, X, Tc, spest, '--')
plt.xlabel('T (ms)')
plt.ylabel('X')
plt.title(f'Speed = {speedest:.2f}')

# Plot some wave profiles
plt.figure(3)
plt.plot(X, sol.y[:N,400], X, sol.y[:N,600])
plt.xlabel('X')
plt.ylabel('V')
plt.legend(['T=20 ms', 'T=30 ms'])


