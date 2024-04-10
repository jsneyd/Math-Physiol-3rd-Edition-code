
#     -------------------------------------------------------------------
# 
#      Adaptation in a hormone receptor, shown by pulsatile stimulation 
#      at different frequencies.
# 
#      For Chapter 16, Section 16.8.1 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def rhs(y, t):
    B, D = y
    H = getH(t)
    dBdt = k1 * H * (1 - B - 2 * D) - km1 * B + 2 * km2 * D - 2 * km2 * B**2
    dDdt = -2 * km2 * D + 2 * km2 * B**2
    return [dBdt, dDdt]

def rhsq(y, t):
    B, D, _ = y
    H = getH(t)
    dBdt = k1 * H * (1 - B - 2 * D) - km1 * B + 2 * km2 * D - 2 * km2 * B**2
    dDdt = -2 * km2 * D + 2 * km2 * B**2
    dQdt = D
    return [dBdt, dDdt, dQdt]

def getH(t):
    if t % T < w:
        return h
    return 0


# Define parameters
k1 = 0.5
k2 = 0.1
km1 = 0.25
km2 = 0.1
stimweight = 0.1

# Compute steady-state response to constant input
T = 1
w = 1
h = stimweight
alph = 2 * k1 * k2 * h / km2
bet = k1 * h + km1
gam = - k1 * h
disc = np.sqrt(bet**2 - 4 * alph * gam)
root1 = (-bet + disc) / (2 * alph)
steady_response = root1**2 * k2 / km2

# Compute response to oscillatory input
n = 40
response = np.zeros(n)
logfreq = np.linspace(-3, 0, n)
freq = 10**logfreq
period = 1 / freq
init = [0, 0]

for i in range(n):
    T = period[i]
    h = 1
    w = stimweight * T / h
    tspan = np.linspace(0, 30 * T, 1000)
    y = odeint(rhs, init, tspan, hmax = w/200)
    init = [y[-1, 0], y[-1, 1]]
    tspan = np.linspace(tspan[-1], tspan[-1] + 2 * T, 1000)
    y = odeint(rhsq, init + [0], tspan, hmax = w/200)
    response[i] = y[-1, 2] / (2 * T)

# Plot results
plt.figure()
plt.semilogx(freq, np.ones(n) * steady_response, '--', label='pulsatile input')
plt.semilogx(freq, response, label='steady input')
plt.xlabel('Frequency of stimulation (1/T)')
plt.ylabel('Average response over period')
plt.legend()



