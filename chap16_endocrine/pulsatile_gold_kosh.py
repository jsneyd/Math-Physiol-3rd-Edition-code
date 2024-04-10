
#     -------------------------------------------------------------------
# 
#      Pulsatile input to the Goldbeter-Koshland equations.
# 
#      For Chapter 16, Section 16.8.2 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def GK(H):
    alpha = 1 - H / v2
    beta = H * (K2 + 1) / v2 - (1 - K2)
    gamma = -K1
    return (-beta + np.sqrt(beta**2 - 4 * alpha * gamma)) / (2 * alpha)

def rhs(y, t):
    s = y[0]
    sst = 1 - s
    H = getH(t)
    dsdt = v2 * sst / (K2 + sst) - H * s / (K1 + s)
    return [dsdt]

def rhsq(y, t):
    s = y[0]
    sst = 1 - s
    H = getH(t)
    dsdt = v2 * sst / (K2 + sst) - H * s / (K1 + s)
    return [dsdt, sst]

def getH(t):
    if t % T < w:
        return h
    return 0

# Define parameters
K1 = 1/100
K2 = 1/100
E2T = 5/100
v2 = E2T
stimweight = 35 / 1000

# Compute steady-state response to constant input
steady_response = 1 - GK(stimweight)

# Compute response to oscillatory input
n = 100
response = np.zeros(n)
logfreq = np.linspace(-2.4, 0, n)
freq = 10**logfreq
period = 1 / freq
init = [0]

for i in range(n):
    T = period[i]
    w = 1
    h = stimweight * T / w
    tspan = np.linspace(0, 50 * T, 1000)  # long enough to sit on the periodic solution
    y = odeint(rhs, init, tspan, hmax = 0.1)
    
    init = [y[-1, 0]]
    tspan = np.linspace(tspan[-1], tspan[-1] + 2 * T, 1000) # integrate over two more periods
    y = odeint(rhsq, init + [0], tspan, hmax = 0.1)
    response[i] = y[-1, 1] / (2 * T)

# Plot results
plt.figure()
plt.semilogx(freq, np.ones(n) * steady_response, '--', label='steady input')
plt.semilogx(freq, response, label='pulsatile input')
plt.xlabel('Frequency (1/T)')
plt.ylabel('Average response ')
plt.legend()

