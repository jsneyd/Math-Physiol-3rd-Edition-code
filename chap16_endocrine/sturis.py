
#     -------------------------------------------------------------------
# 
#      Solve the Sturis model for pulsatile insulin secretion.
# 
#      For Chapter 16, Section 16.7.3 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Set parameters
input_ = 216
Vp = 3
Vi = 11
Vg = 10
E = 0.2
tp = 6
ti = 100
td = 12
Rm = 209
a1 = 6.6
C1 = 300
C2 = 144
C3 = 100
C4 = 80
C5 = 25.86
Ub = 72
U0 = 4
Um = 94
Rg = 180
alpha = 7.5
beta = 1.77


initial = [1, 1, 1, 1, 1, 1]
tspan = np.linspace(0, 2000, 1000)

# Function to solve the ODE system
def sturisfun(y, t):
    Ip, Ii, G, h1, h2, h3 = y
    
    # Define the ODEs
    dydt = [
        eff1(G) - (Ip / Vp - Ii / Vi) * E - Ip / tp,
        (Ip / Vp - Ii / Vi) * E - Ii / ti,
        eff5(h3) + input_ - eff2(G) - eff3(G) * eff4(Ii),
        (Ip - h1) / td,
        (h1 - h2) / td,
        (h2 - h3) / td
    ]
    
    return dydt

# Define efficiency functions
def eff1(z):
    return Rm / (1 + np.exp(-z / (C1 * Vg) + a1))

def eff2(z):
    return Ub * (1 - np.exp(-z / (C2 * Vg)))

def eff3(z):
    return z / (C3 * Vg)

def eff4(z):
    kappa = (1 / Vi + 1 / (E * ti)) / C4
    return U0 + (Um - U0) / (1 + (kappa * z) ** (-beta))

def eff5(z):
    return Rg / (1 + np.exp(alpha * (z / (C5 * Vp) - 1)))

# Solve the ODE system
sol = odeint(sturisfun, initial, tspan)

# Plot results
plt.figure(1)
plt.plot(tspan, sol[:, 0] / Vp)
plt.xlabel('Time (minutes)')
plt.ylabel('Insulin (mU/ml)')
plt.figure(2)
plt.plot(tspan, sol[:, 2] / (10 * Vg))  # Output per dL, as in Sturis
plt.xlabel('Time (minutes)')
plt.ylabel('Glucose (mg/dl)')
plt.show()
