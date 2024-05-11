
# -------------------------------------------------------------------
# 
#  Simulate cochlear waves with the Lesser-Berkley model
#  
#  For Chapter 20, Section 20.2.3 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#  Written by James Keener and James Sneyd
# 
# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# Set parameters
num = 200
L = 3.5  # Units of cm
l = 0.35
lam = 1.5
xi = np.linspace(0, 1, num)
k = (1.e7) * np.exp(-lam * xi * L)  # don't forget that x=xi*L
mass = 0.05
r = 3000 * np.exp(-lam * xi * L)
# Two cases
w_list = [800, 1500]

for w in w_list:
    Z = 1j * w * mass + r - 1j * k / w
    W = -1j * Z / (w * L)
    sigma = l / L

    # Initialize
    N = 100
    alpha = np.zeros((N + 1, N + 1), dtype=np.complex_)
    f = np.zeros(N + 1, dtype=np.complex_)
    A = np.zeros(N + 1, dtype=np.complex_)

    # Get matrix coefficients
    for n in range(N + 1):
        for m in range(N + 1):
            alpha[m, n] = 2 * np.cosh(n * np.pi * sigma) * np.trapz(np.cos(n * np.pi * xi) * np.cos(m * np.pi * xi) / W, xi)
            f[m] = -np.trapz(xi * (2 - xi) * np.cos(m * np.pi * xi) / W, xi)

        alpha[n, n] = alpha[n, n] + n * np.pi * np.sinh(n * np.pi * sigma) / 2

    f[0] = f[0] - sigma

    # Solve linear equations for Fourier coefficients
    A = np.linalg.solve(alpha, f)

    # Construct phi on y=0
    y = 0.0  # Evaluate eta on the membrane
    phi = xi * (1 - xi / 2) - sigma * y * (1 - y / (2 * sigma))
    for n in range(N + 1):
        phi = phi + A[n] * np.cosh(n * np.pi * (sigma - y)) * np.cos(n * np.pi * xi)

    Fhat = 1  # Driving force
    phi = phi / (1j * w * L * Fhat)
    eta = 2 * phi / W

    # Plot results
    plt.figure(2 * w_list.index(w) + 1)
    plt.plot(xi * L, np.real(eta), 'r')
    plt.plot(xi * L, np.abs(eta), '--b')
    plt.plot(xi * L, -np.abs(eta), '--b')
    plt.xlabel('x (cm)')
    plt.ylabel('amplitude')
    plt.xlim([0, L])
    plt.title(f'$\omega = {w}$/s')

    # Animate the wave
    N = 200
    times = np.linspace(0, 10 / w, N)
    plt.figure(2 * w_list.index(w) + 2)
    for j in range(N):
        etawave = eta * np.exp(1j * w * times[j])
        plt.plot(xi * L, np.real(etawave), 'r')
        plt.plot(xi * L, np.abs(etawave), '--b')
        plt.plot(xi * L, -np.abs(etawave), '--b')
        plt.xlabel('x (cm)')
        plt.ylabel('amplitude')
        plt.title(f'$\omega = {w}$/s')
        plt.xlim([0, L])
        plt.draw()
        plt.pause(0.02)

plt.show()
