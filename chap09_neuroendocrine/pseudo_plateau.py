#  ------------------------------------
#  Pseudo-plateau bursting. A modification of the Chay-Keizer model.
#  Modified from the original code of Teka, Bertram, etc.

#  Used for Keener and Sneyd, Mathematical Physiology, Chapter 9, Section 9.4.

#  Variables:
#     V -- membrane potential
#     n -- delayed rectifier activation variable
#     c -- cytosolic calcium concentration
#  ------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# ODE equations
def deRHS(sol, t):
    global gkatp, vk, vca, cm, gk, gca, vm, sm, vn, sn, taun, gkca, kd, kpmca, f, alpha

    # There are three variables
    v, n, c = sol

    # Steady state functions
    ninf = 1 / (1 + np.exp((vn - v) / sn))
    minf = 1 / (1 + np.exp((vm - v) / sm))

    # Ikca
    Ikca = gkca / (1 + (kd / c)**3) * (v - vk)

    # Calcium Handling

    # ICa
    Ica = gca * minf * (v - vca)

    # Ik
    Ik = gk * n * (v - vk)

    # Ikatp
    Ikatp = gkatp * (v - vk)

    # Ca fluxes
    Jmem = -(alpha * Ica + kpmca * c)

    # Equations
    dv = -(Ik + Ica + Ikca + Ikatp) / cm
    dn = (ninf - n) / taun
    dc = f * Jmem

    return [dv, dn, dc]

# Global parameters
gkatp = 180
vk = -75
vca = 25
cm = 5300
gk = 2700
gca = 1000
vm = -20
sm = 12
vn = -12
sn = 5
taun = 18.7
gkca = 400
kd = 0.3
kpmca = 0.18
f = 0.01
alpha = 4.50e-6

# Initial conditions
v0 = -65
n0 = 0
c0 = 0.1
init = [v0, n0, c0]

total = 10000
tstep = 0.01
# Specify the output points
tspan = np.arange(0, total + tstep, tstep)

# Solve the ODE
sol = odeint(deRHS, init, tspan)

plt.figure(1)
plt.plot(tspan / 1000, sol[:, 0])
plt.ylabel('V (mV)')
plt.xlabel('t (s)')

plt.figure(2)
plt.plot(sol[:, 2], sol[:, 0])
plt.ylabel('V (mV)')
plt.xlabel('Ca$^{2+}$')

plt.show()



