
# -------------------------------------------------------------------
# 
#  The fire-diffuse-fire model.
# 
#  For Chapter 7, Section 7.8.1 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#  Written by James Keener and James Sneyd
# 
# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

betalist = np.array([0.5, 0.1, 0.01])
dkappa = 0.001
eta = np.arange(0.01, 50, 0.05)
J = len(betalist)

G = np.zeros((J, len(eta)))
K = np.zeros(J)

for j in range(J):
    beta = betalist[j]
    g = np.zeros((1000, len(eta)))
    
    for n in range(1, 1001):
        g[n-1, :] = np.sqrt(1 / (4 * np.pi * n * eta)) * np.exp(-n / (4 * eta) - beta**2 * n)
    
    G[j, :] = np.sum(g, axis=0)
    K[j] = beta

plt.figure(1)
for j in range(J):
    plt.plot(eta, G[j, :], linewidth=2, label=f'beta={betalist[j]}')
plt.xlabel(r'$\eta$', fontsize=16)
plt.ylabel(r'$g_\beta(\eta)$', fontsize=16)
plt.legend()

plt.figure(2)
for j in range(J):
    plt.plot(G[j, :], eta, linewidth=2, label=f'beta={betalist[j]}')
plt.ylabel('delay', fontsize=16)
plt.xlabel('threshold', fontsize=16)
plt.axis([0, 1, 0, 5])
plt.legend()

plt.figure(3)
plt.semilogy(np.sqrt(betalist), np.max(G, axis=1), linewidth=2, label='max(G)')
plt.semilogy(np.sqrt(betalist), np.exp(-np.sqrt(betalist)), '--', linewidth=2, label='exp(-sqrt(beta))')
plt.xlabel(r'$\beta$', fontsize=16)
plt.ylabel(r'$g_{max}$', fontsize=16)
plt.legend()


