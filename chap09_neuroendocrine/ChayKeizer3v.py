#   -------------------------------------------------------------------
# 
#    Code for the three-variable Chay-Keizer model of electrical bursting 
#    in pancreatic beta cells.
# 
#    For Chapter 9, Section 9.1.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

from numpy import exp, arange
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def deRHS(u, t):
    v, n, ca = u
    # Gating functions
    alphamca = -0.1 * ((v + vprime) - 25) / (exp(-((v + vprime) - 25) / 10) - 1)
    betamca = 4 * exp(-(v + vprime) / 18)
    alphahca = 0.07 * exp(-(v + vprime) / 20)
    betahca = 1 / (exp(-((v + vprime) - 30) / 10) + 1)
    alphan = -0.01 * ((v + vstar) - 10) / (exp(-((v + vstar) - 10) / 10) - 1)
    betan = 0.125 * exp(-(v + vstar) / 80)

    minf = alphamca / (alphamca + betamca)
    hinf = alphahca / (alphahca + betahca)

    # Ionic currents
    ica = gca * minf**3 * hinf * (v - vca)
    ik = gk * n**4 * (v - vk)
    ikca = gkca * ca / (ca + Kd) * (v - vk)
    il = gl * (v - vl)

    # Differential equations
    vp = -1 / cm * (ica + ik + ikca + il)
    np = lam * (alphan * (1 - n) - betan * n)
    cap = f * (-alpha * ica - kca * ca)

    return [vp, np, cap]


# Global parameters
lam, vca, vk, gk, gca, vprime, vstar, gkca, Kd, vl, cm, f, alpha, kca, gl = [
    0.5, 100, -75, 12, 13, 50, 30, 0.09, 1, -40, 4.5, 0.01, 0.01525, 0.04, 0.04
]

# Simulation parameters
total = 10000
tstep = 1

# Specify the output points
tspan = arange(0, total + tstep, tstep)

# Initial conditions
v0, n0, ca0 = -54.77, 0.00044, 0.1075
u0 = [v0, n0, ca0]

# Solve the ODE system
sol = odeint(deRHS, u0, tspan)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(sol[:, 2], sol[:, 0])
plt.xlabel('Ca++')
plt.ylabel('V')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(tspan / 1000, sol[:, 0])
plt.xlabel('t (s)')
plt.ylabel('V')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(tspan / 1000, sol[:, 2])
plt.xlabel('t (s)')
plt.ylabel('Ca++')
plt.show()
