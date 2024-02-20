
#   -------------------------------------------------------------------
# 
#    ODE integrator for the MNT Purkinje fiber model equations.
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
def MNTdeRHS(s, t):
    V = s[0]
    alpbet = ab(V)
    currents = IV(s)
    Fv = -(sum(currents) + Inp) / cm
    Fg = alpbet[::2] * (1 - s[1:]) - alpbet[1::2] * s[1:]
    return np.concatenate(([Fv], Fg))


def IV(S):
    V = S[0]
    m = S[1]
    h = S[2]
    d = S[3]
    f = S[4]
    q = S[5]
    r = S[6]
    s = S[7]
    x1 = S[8]
    x2 = S[9]
    
    gna = gnabar * m ** 3 * h
    INa = gna * (V - Vna)
    dpr = 1 / (1 + np.exp(-0.15 * (V + 40)))
    Isi = (0.8 * d * f + 0.04 * dpr) * (V - Vsi)
    IK2bar = (np.exp(0.04 * (V + 110)) - 1) / (np.exp(0.08 * (V + 60)) + np.exp(0.04 * (V + 60)))
    IK2 = 2.8 * IK2bar * s
    Ix1 = 1.2 * x1 * (np.exp(0.04 * (V + 95)) - 1) / np.exp(0.04 * (V + 45))
    Ix2 = x2 * (25 + 0.385 * V)
    Icl = 2.5 * q * r * (V - Vcl)
    IK1 = IK2bar + 0.2 * (V + 30) / (1 - np.exp(-0.04 * (V + 30)))
    Inab = 0.105 * (V - 40)
    Iclb = 0.01 * (V + 70)

    return np.array([INa, Isi, IK2, Ix1, Ix2, Icl, IK1, Inab, Iclb])



def ab(V):  
    # there are 5 parameters C for each alpha, beta
    alphabeta = (C[:,0]*np.exp(C[:,1]*(V-V0))+C[:,2]*(V-V0))/(1+C[:,3]*np.exp(C[:,4]*(V-V0)))
    return alphabeta


# Parameters
gnabar = 150.
Vna = 40.
Vcl = -70
Vsi = 70

cm = 10.

# The matrix C has coefficients for alpha beta
C = np.array([[0, 0, 1, -1, -0.1],
              [40, -0.056, 0, 0, 0],
              [0.0085, -0.184, 0, 0, 0],
              [2.5, 0, 0, 1, -0.082],
              [0, 0, 0.002, -1, -0.1],
              [0.02, -0.089, 0, 0, 0],
              [0.000987, -0.04, 0, 0, 0],
              [1, 0, 0, 1, -0.087],
              [0, 0, 0.008, -1, -0.1],
              [0.08, -0.0888, 0, 0, 0],
              [0.00018, -0.04, 0, 0, 0],
              [0.02, 0, 0, 1, -0.087],
              [0, 0, 0.001, -1, -0.2],
              [5.e-5, -0.067, 0, 0, 0],
              [0.0005, 0.083, 0, 1, 0.057],
              [0.0013, -0.06, 0, 1, -0.04],
              [0.083e-4, 0, 0, 1, -0.2],
              [0.0003, -0.06, 0, 1, -0.04]])

V0 = np.array([-47, -72, -71, -10, -40, -40, -60, -26, 0, 0, -80, -26, -52, -52, -50, -20, -19, -20])

tstep = 0.1
t_end = 5000

# Initial conditions
Inp = 0
V = -50
s0 = [V, 0.5, 0, 0.21, 0, 0, 0, 1, 0.15, 0]

# Time points
tspan = np.arange(0, t_end, tstep)


# Calculate the solution
S = odeint(MNTdeRHS, s0, tspan)

# Plotting
T = tspan
V = S[:, 0]
m = S[:, 1]
h = S[:, 2]
d = S[:, 3]
f = S[:, 4]
q = S[:, 5]
r = S[:, 6]
s = S[:, 7]
x1 = S[:, 8]
x2 = S[:, 9]

plt.figure(1)
plt.subplot(2, 1, 1)
plt.plot(T, V, linewidth=2)
plt.xlabel('t (ms)', fontsize=20)
plt.ylabel('V (mV)', fontsize=20)

plt.subplot(2, 1, 2)
plt.plot(T, m, 'r', T, h, 'b', T, d, 'g', linewidth=2)
plt.legend(['m', 'h', 'd'], fontsize=16)
plt.xlabel('t (ms)', fontsize=20)

plt.figure(2)
plt.plot(T, f, T, q, T, r, T, s, T, x1, T, x2)
plt.legend(['f', 'q', 'r', 's', 'x1', 'x2'], fontsize=16)

INa=np.zeros(len(T))
Isi=np.zeros(len(T))
IK2=np.zeros(len(T))
Ix1=np.zeros(len(T))
Ix2=np.zeros(len(T))
INx2=np.zeros(len(T))
Icl=np.zeros(len(T))
IK1=np.zeros(len(T))
Inab=np.zeros(len(T))
Iclb=np.zeros(len(T))

for j in range(len(T)):
    s = S[j,:]
    currents = IV(s)
    INa[j] = currents[0]
    Isi[j] = currents[1]
    IK2[j] = currents[2]
    Ix1[j] = currents[3]
    Ix2[j] = currents[4]
    Icl[j] = currents[5]
    IK1[j] = currents[6]
    Inab[j] = currents[7]
    Iclb[j] = currents[8]

plt.figure(3)
plt.plot(T, Isi, T, IK2, T, Ix1, T, Ix2, T, IK1, T, Inab, T, Iclb, linewidth=2)
plt.legend(['Isi', 'IK2', 'Ix1', 'Ix2', 'IK1', 'Inab', 'Iclb'], fontsize=16)
plt.xlabel('t (ms)', fontsize=20)

plt.figure(4)
plt.plot(T, INa)
plt.legend(['Ina'], fontsize=16)
plt.xlabel('t (ms)', fontsize=20)
plt.show()
