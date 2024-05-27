#   -------------------------------------------------------------------
# 
#    Code to calculate and plot an effective diffusion coefficient.
# 
#    For Chapter 8, Figure 8.21 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

def deff(l, delta):
    N = 50
    A = np.zeros((N, N))
    R = np.zeros(N)
    P = np.zeros(N)

    for kj in range(1, N + 1):
        for n in range(1, N + 1):
            k = 1 + 2 * (kj - 1)
            A[kj - 1, n - 1] = (2 * n * np.tanh(2 * n * np.pi * l * delta) / (k * np.tanh(k * np.pi * l * (1 - delta))) + 1) * n / (4 * n**2 - k**2)
        R[kj - 1] = 1 / (2 * np.pi * k**2)
        P[kj - 1] = np.tanh(2 * kj * np.pi * l * delta)
    C = np.linalg.solve(A, R)
    return np.dot(P, C) / l + delta

# Plotting the first set of data
dellist = [0.2, 0.1, 0.01]
llist = np.arange(0.01, 4.01, 0.01)
for delta in dellist:
    Deff = np.array([deff(l, delta) for l in llist])
    plt.figure(1)
    plt.plot(llist, Deff)
    plt.xlabel('l/L')
    plt.ylabel('D_{eff}/D')
    plt.legend([f'\u0394={delta}' for delta in dellist])

plt.text(0.8, 0.5, '\u0394=0.2', fontsize=18)
plt.text(1.2, 0.3, '\u0394=0.1', fontsize=18)
plt.text(2, 0.135, '\u0394=0.01', fontsize=18)
plt.show()

# Plotting the second set of data
dellist = np.arange(0.001, 1.01, 0.01)
llist = [0.00001, 0.2, 0.5]
for lidx, l in enumerate(llist):    
    Deff = np.array([deff(l, delta) for delta in dellist])
    plt.figure(2)
    plt.plot(dellist, Deff)
    plt.xlabel('\u0394')
    plt.ylabel('D_{eff}/D')
    plt.legend([f'l/L={l}' for l in llist])

    if lidx == 0:
        p = dellist * (Deff - 1)
        q = np.polyfit(dellist[5:-1], p[5:-1], 1)
        pfit = np.polyval(q, dellist)

        plt.figure(3)
        plt.plot(dellist, p, dellist, pfit, '--')

plt.text(0.1, 0.97, 'l/L=0', fontsize=18)
plt.text(0.08, 0.8, 'l/L=0.2', fontsize=18)
plt.text(0.18, 0.5, 'l/L=0.5', fontsize=18)
plt.show()
