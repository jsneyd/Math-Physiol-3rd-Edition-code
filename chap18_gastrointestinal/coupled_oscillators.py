
#     -------------------------------------------------------------------
# 
#      Simulate phase equations.
# 
#      For Chapter 18, Section 18.5.2 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# parameters
N = 32  # number of coupled oscillators
x1mG1 = 20
Omega = -10 / 31

# natural frequencies
num = np.arange(1, N + 1)
nat_freq = x1mG1 + Omega * (num - 1)

# pick a value of delta
deltalist = [32, 18]
for delta in deltalist:
    # coupling matrix
    K = -2 * np.eye(N) + np.diag(np.ones(N - 1), -1) + np.diag(np.ones(N - 1), 1)

    # set up the differential equation solve
    tspan = (0, 2000)  # time span
    t_eval = np.arange(tspan[0], tspan[1] + 1, 1)  # time points for evaluation

    # initial data for integration
    s0 = np.zeros(N + 1)

    # define the differential equation
    def deRHS(t, s):
        f = np.sin(s[1:])
        F = Omega + delta * K.dot(f)
        F1 = x1mG1 + delta * f[0]
        return np.concatenate(([F1], F))

    # integrate the differential equation. odeint is no good because 
    # you need to use a stiff method.
    sol = solve_ivp(deRHS, tspan, s0, method='Radau', t_eval=t_eval)

    plt.figure(3 * deltalist.index(delta) + 2)
    plt.plot(sol.y[1:, -1] / tspan[-1], '*', linewidth=2)
    plt.xlabel('Oscillator #', fontsize=20)
    plt.ylabel('Average Phase Difference')
    plt.title(f'$\delta = {delta}$', fontsize=18)
    
    plt.figure(3 * deltalist.index(delta) + 3)
    Cs = np.cumsum(sol.y[:, -1]) / tspan[-1]
    plt.plot(num, Cs[1:], '*', num, nat_freq, '--')
    plt.xlabel('Oscillator #', fontsize=20)
    plt.ylabel('Frequency')
    plt.title(f'$\delta = {delta}$', fontsize=18)


