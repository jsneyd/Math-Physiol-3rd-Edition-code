#   -------------------------------------------------------------------
# 
#    Code for the five-variable Chay-Keizer model of electrical bursting 
#    in pancreatic beta cells.
# 
#    For Chapter 9, Section 9.1.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

from numpy import exp,arange
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# ODE system for the Chay-Keizer model
def deRHS(u, t):
    v, mca, hca, n, ca = u

    # Gating functions
    alphamca = -0.1 * ((v + vprime) - 25) / (exp(-((v + vprime) - 25) / 10) - 1)
    betamca = 4 * exp(-(v + vprime) / 18)
    alphahca = 0.07 * exp(-(v + vprime) / 20)
    betahca = 1 / (exp(-((v + vprime) - 30) / 10) + 1)
    alphan = -0.01 * ((v + vstar) - 10) / (exp(-((v + vstar) - 10) / 10) - 1)
    betan = 0.125 * exp(-(v + vstar) / 80)

    # Ionic currents
    ica = gca * mca**3 * hca * (v - vca)
    ik = gk * n**4 * (v - vk)
    ikca = gkca * ca / (ca + Kd) * (v - vk)
    il = gl * (v - vl)

    # Differential equations
    vp = -1 / cm * (ica + ik + ikca + il)
    mcap = alphamca * (1 - mca) - betamca * mca
    hcap = alphahca * (1 - hca) - betahca * hca
    np = alphan * (1 - n) - betan * n
    cap = f * (-k1 * ica - kca * ca)

    return [vp, mcap, hcap, np, cap]

# Global parameters
vprime, vstar, gca, vca, gk, gkca, Kd, vk, gl, vl, cm, k1, kca, f = [
    50, 30, 13, 100, 12, 0.09, 1, -75, 0.04, -40, 4.5, 0.01525, 0.04, 0.01
]

# Simulation parameters
total = 15000
tstep = 1

# Specify the output points
tspan = arange(0, total + tstep, tstep)

# Initial conditions
v0, mca0, hca0, n0, ca0 = -54.774, 0.027532, 0.086321, 0.00044035, 0.10749
u0 = [v0, mca0, hca0, n0, ca0]

# Solve the ODE system
sol = odeint(deRHS, u0, tspan)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(sol[:, 4], sol[:, 0])
plt.xlabel('Ca++')
plt.ylabel('V')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(tspan / 1000, sol[:, 0])
plt.xlabel('t (s)')
plt.ylabel('V')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(tspan / 1000, sol[:, 4])
plt.xlabel('t (s)')
plt.ylabel('Ca++')
plt.show()
