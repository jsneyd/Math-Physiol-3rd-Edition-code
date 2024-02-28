
#   -------------------------------------------------------------------
# 
#    Code to examine phase synchrony
# 
#    For Chapter 12, Section 12.4.4 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# Set parameters
N = 30  # Size of grid
Xcent = 12
Ycent = 18  # Center of the grid
scale = 50
amp = 0.3

# Create a grid of random numbers
X, Y = np.meshgrid(np.arange(1, N+1), np.arange(1, N+1))
rad = (X-Xcent)**2 + (Y-Ycent)**2
P = np.random.rand(N, N)
Pf = 2 * P * (1 - P)

smth = np.exp(-rad/scale) + Pf

# Normalize the grid
smth = smth / np.max(smth)

# Plot the original grid
plt.figure(1)
plt.pcolor(smth)

# Create the coupling matrix
offdiag = np.ones(N-1)
for j in range(N-1):
    offdiag = np.concatenate([offdiag, [0], np.ones(N-1)])


A = np.diag(offdiag, 1) + np.diag(offdiag, -1) + np.diag(np.ones(N*(N-1)), N) + np.diag(np.ones(N*(N-1)), -N)
Ad = np.sum(A, axis=1)
A = A - np.diag(Ad)

# Create the right-hand side
rhs = smth.flatten()

# Project out the adjoint nullspace
srhs = np.sum(rhs)
rhs = rhs - np.ones(N**2) * srhs / N**2

# Replace last row of A by ones
A[N**2-1] = np.ones(N**2)
rhs[N**2-1] = 1

# Solve the linear system of equations
phase = np.linalg.solve(A, rhs)
nphase = -phase.reshape(N, N)

# Plot the phase distribution
plt.figure(2)
plt.pcolor(nphase)

