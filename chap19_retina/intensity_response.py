
#     -------------------------------------------------------------------
# 
#      Intensity-response curves in the model of light adaptation in
#      in cones.
# 
#      For Chapter 19, Section 19.2.2 and Exercise 19.3 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------


#  The IR curves are computed by first holding the cone at a fixed
#  background light level (of I0_base) so that the cone reaches steady state
#  for this background. Then the background is either increased or decreased
#  and the peak response (either the max or the min) is recorded.

#  Because we are doing the perturbation using the background light levels,
#  stim is always held at zero. This avoids complications due to the kernel
#  getting too large, and the light stimulus going negative. Which would be
#  bad.


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def rhs(t, U, par, stim):
    p, x, y, z, v = U
    kernel = (par['eta'] / par['tau1'] / 6) * ((t / par['tau1'])**3) * np.exp(-t / par['tau1'])
    s = par['eta'] * par['I0'] + stim * kernel

    phi = getphi(y, par)

    dUdt = [0, 0, 0, 0, 0]
    dUdt[0] = s * (1 - p) - par['k1'] * p
    dUdt[1] = phi - (par['gam'] - par['delta']) * x * p - par['delta'] * x
    dUdt[2] = ((x**3) * np.exp(-v) - y) / par['tauy']
    dUdt[3] = (((1 - par['kappa']) / (1 + par['kappa'])) * (x**3) * np.exp(-v)
               + (2 * par['kappa'] / (1 + par['kappa'])) * y - z) / par['tauz']
    dUdt[4] = ((x**3) * np.exp(-v) - ((1 + par['kappa']) / 3) * z
               + (par['kappa'] / 2) * y + ((4 + par['kappa']) / (6 * par['vK'])) * (v - par['vK'])) / par['taum']
    
    return dUdt

def getphi(y, par):
    v = par['vK'] * (1 - y)
    x = (y * np.exp(v))**(1 / 3)
    I = (np.exp(-v / par['s1']) - 1) / par['s2']
    p = par['eta'] * I / (par['k1'] + par['eta'] * I)

    # Comment out the one you don't want
    phi = x * (par['delta'] + (par['gam'] - par['delta']) * p)
    # phi = 4 + 84 / (1 + (y / 0.34) ** 4)  # Using A(y) instead 
    return phi


par = {'Vstar': 35.7, 's1': 1.59 / 35.7, 's2': 1130, 'vK': -13 / 35.7, 'tau1': 0.012,
       'taum': 0.02, 'tauz': 0.04, 'tauy': 0.07, 'k1': 35.4, 'gam': 303, 'delta': 5,
       'kappa': 0.1, 'eta': 52.5, 'I0': 0}

baselist = [0.0001, 0.001, 0.01, 0.1, 1]
I0log = np.linspace(-5, 1, 100)
vkeep = I0log.reshape(-1, 1)

for j in range(5):
    I0_base = baselist[j]

    I0list = 10**I0log
    res = []

    # First find the steady state
    IC = [0, 1, 1, 1, 0]
    tspan = np.linspace(0, 5, 200)
    par['I0'] = I0_base
    sol = solve_ivp(lambda t, U: rhs(t, U, par, 0), [0, 5], IC, t_eval=tspan)
    IC = sol.y[:, -1]  # initial condition for all subsequent runs

    tspan = np.linspace(0, 2, 200)
    for I0 in I0list:
        par['I0'] = I0
        sol = solve_ivp(lambda t, U: rhs(t, U, par, 0), [0, 2], IC, t_eval=tspan)
        U = sol.y
        if I0_base < I0:
            res.append(np.min(U[4]))
        else:
            res.append(np.max(U[4]))

    res = np.array(res)
    plt.plot(I0log, -res * par['Vstar'])
    plt.plot(np.log10(I0_base), -par['Vstar'] * IC[4], 'ro')  # plot the steady-state value for I0_base
    plt.xlabel('log(I_0)')
    plt.ylabel('-V')

    # Add the Naka-Rushton equation
    nakarushton = 16 * I0list / (I0list + 0.015)
    plt.plot(I0log, nakarushton, '--k', linewidth=3)

plt.show()

