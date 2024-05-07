
#     -------------------------------------------------------------------
# 
#      Model of light adaptation in
#      in cones.
# 
#      For Chapter 19, Section 19.2.2 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def rhs(t, X):
    p, x, y, z, v = X
    kernel = (eta/tau1/6)*((t/tau1)**3)*np.exp(-t/tau1)
    s = eta*I0 + stim*kernel
    phi = (y * np.exp(vK * (1 - y))) ** (1 / 3) * (delta + (gam - delta) * eta * (np.exp(-vK * (1 - y) /
        s1) - 1) / (s2 * k1 + eta * (np.exp(-vK * (1 - y) / s1) - 1)))

    ppdot = s * (1 - p) - k1 * p
    xdot = phi - (gam - delta) * x * p - delta * x
    ydot = ((x**3) * np.exp(-v) - y) / tauy
    zdot = (((1 - kappa) / (1 + kappa)) * (x**3) * np.exp(-v) + (2 * kappa / (1 + kappa)) * y - z) / tauz
    Vdot = ((x**3) * np.exp(-v) - ((1 + kappa) / 3) * z + (kappa / 2) * y + ((4 + kappa) / (6 * vK)) * (v - vK)) / taum

    return [ppdot, xdot, ydot, zdot, Vdot]

vstar = 35.7
s1 = 1.59 / vstar
s2 = 1130
vK = -13 / vstar
tauy = 0.07
k1 = 35.4
gam = 303
delta = 5
kappa = 0.1
eta = 52.5
taum = 0.016
tauz = 0.04
tau1 = 0.012
eta = 52.5

stimlist = [0, 0.0001, 0.001, 0.01, 0.1]
IC = [0, 1, 1, 1, 0]
I0 = 0

for j in range(5):
    stim = stimlist[j]
    tspan = np.linspace(0, 2, 200)
    sol = solve_ivp(rhs, [0, 2], IC, t_eval=tspan)
    t = sol.t + 2 * j
    plt.plot(t, sol.y[4,:] * vstar,'r')
    IC = sol.y[:, -1]
