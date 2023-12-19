#   -------------------------------------------------------------------
# 
#    This computes the fast phase plane of the Hodgkin-Huxley model,
#    with n and h fixed.
# 
#    For Chapter 5, Fig. 5.7 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def deRHS(s, t):
    V = s[0] - Veq
    m = s[1]
    h = s[2]
    n = s[3]

    Input = Inp / np.cosh(np.mod(t, BCL) - BCL / 10) + Iapp

    AM = 0.1 * (25. - V) / (np.exp(0.1 * (25. - V)) - 1.)
    BM = 4. * np.exp(-V / 18)
    Fm = AM * (1. - m) - BM * m
    Fh = 0
    Fn = 0
    gna = gnabar * m ** 3 * h
    gk = gkbar * n ** 4
    Icl = gl * (V - Vl)

    Fv = -gna * (V - Vna) - gk * (V - Vk) - Icl
    return [Fv + Input, Fm, Fh, Fn]
    
    
# Set default plot parameters
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['lines.linewidth'] = 2.0

format_spec_f = '%6.2f\n'


# Define global parameters
Veq = 0  # Original HH formulation
# p.Veq = -65  # Physiological formulation
Tfact = 1  # Correction factor for temperatures other than 6.3°C
# Tfact = 0.977  # Corresponds to 0°C
# Tfact = 1.085  # Corresponds to 30°C

gnabar = 120.0
gkbar = 36.0
gl = 0.3
# p.gkbar=0  # to shut down currents
# p.gl = 0
Vna0 = 115
Vk0 = -12
Vl0 = 10.5988

Vna = Tfact * (Vna0 + Veq) - Veq
Vk = Tfact * (Vk0 + Veq) - Veq
Vl = Tfact * (Vl0 + Veq) - Veq

h0 = 0.5961
n0 = 0.3177
Iapp = -2  # p.Iapp is the amplitude of the steady current stimulus

# Start by plotting nullclines
Vm = np.arange(-20, 120, 0.1)
AM = 0.1 * (25. - Vm) / (np.exp(0.1 * (25. - Vm)) - 1.)
BM = 4. * np.exp(-Vm / 18)

Minf = AM / (AM + BM)
M = np.arange(0, 1, 0.001)
gk = gkbar * n0 ** 4
gna = gnabar * M ** 3 * h0

Vs = (gna * Vna + gk * Vk + gl * Vl + Iapp) / (gna + gk + gl)

# This is the dv/dt=0 curve
plt.figure(1)
plt.plot(Vm, Minf, '--', label='Nullcline m_inf(v)')
plt.plot(Vs, M, '--', label='Nullcline dv/dt=0')
plt.xlabel('V')
plt.ylabel('m')
plt.legend()
plt.title(f'I_app = {Iapp:.2f}')
plt.axis([-10, 120, 0, 1])

# closeup view of the saddle point
plt.figure(2)
plt.plot(Vm, Minf, '--', label='Nullcline m_inf(v)')
plt.plot(Vs, M, '--', label='Nullcline dv/dt=0')
plt.xlabel('V')
plt.ylabel('m')
plt.legend()
plt.title(f'I_app = {Iapp:.2f}')
plt.axis([-10, 10, 0, 0.2])

Inp = 0  # This is the amplitude of the initial input stimulus
# Inp = 3.5
# Inp = 13.4

BCL = 39
# BCL = 0

tstep = 0.01
t_end = 10
V = 6.5098  # The voltage threshold

# Add some trajectories to both graphs
vlist = [10, 7, -5, -5, 8.29099503, -5]
mlist = [0, 0, 0.22, 0.24, 0, 0.2308685]

for j in range(len(vlist)):
    # Initial data
    # Specify the output points
    tspan = np.arange(0, t_end + tstep, tstep)

    s0 = [vlist[j], mlist[j], h0, n0]
    s = odeint(deRHS, s0, tspan)
    V = s[:, 0]
    m = s[:, 1]

    plt.figure(1)
    plt.plot(V, m, linewidth=2)
    plt.figure(2)
    plt.plot(V, m, linewidth=2)





