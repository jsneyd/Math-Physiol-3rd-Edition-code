
#     -------------------------------------------------------------------
# 
#      Hamer model of light adaptation in rods.
# 
#      For Chapter 19, Section 19.3 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def rods(t, X):
    R, E, g, c, b = X

    light = stim * ((t > ton) & (t < toff))
    f = (g / gdark)**ng

    dRdt = light - R / tauR
    dEdt = v * R - E / tauE
    dgdt = Amax / (1 + (c / Kc)**nc) - (beta_dark + betaE * E) * g
    dcdt = alpha * f * Jdark - gamma * (c - cmin) - kon * (bt - b) * c + koff * b
    dbdt = kon * (bt - b) * c - koff * b

    return [dRdt, dEdt, dgdt, dcdt, dbdt]


ng = 3
Jdark = 72
tauR = 0.4
tauE = 2
v = 1160
bt = 660
gamma = 73
kon = 0.1
koff = 0.8
cmin = 0.005
Kc = 0.9
nc = 3
Amax = 10
betaE = 170
beta_dark = 0.13
alpha = 0.8
gdark = 2

IC = [0, 0, 1.1068, 0.7106, 19.3133]
tspan = np.linspace(0, 30, 1000)
options = {'max_step': 0.01}

# First find the dark state, to use as initial condition
ton = 0
toff = 0
stim = 0
sol = solve_ivp(rods, [0, 30], IC, t_eval=tspan, method='Radau', **options)
IC = sol.y[:, -1]

# compute flash responses
ton = 1
toff = 1.01
stimlist = [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000]
for stim in stimlist:
    sol = solve_ivp(rods, [0, 30], IC, t_eval=tspan, method='Radau', **options)
    photocurrent = (sol.y[2] / IC[2])**ng * Jdark
    plt.plot(sol.t, photocurrent[0] - photocurrent)
plt.xlabel('Time (s)')
plt.ylabel('Photocurrent (pA)')

# compute step responses
plt.figure()
tspan = np.linspace(0, 100, 1000)
ton = 1
toff = 60
stimlist = [1e-6, 2e-6, 1e-5, 2e-5, 1e-4, 1e-3, 1e-2]
for stim in stimlist:
    sol = solve_ivp(rods, [0, 100], IC, t_eval=tspan, method='Radau', **options)
    photocurrent = (sol.y[2] / IC[2])**ng * Jdark
    plt.plot(sol.t, photocurrent[0] - photocurrent)
plt.xlabel('Time (s)')
plt.ylabel('Photocurrent (pA)')



