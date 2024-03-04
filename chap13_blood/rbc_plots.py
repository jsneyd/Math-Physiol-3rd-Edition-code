#   -------------------------------------------------------------------
# 
#  Code to solve red blood cell delay differential equation
#  using method of lines.
# 
#    For Chapter 13, Section 13.1.3 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# RHS for time-dependent simulation
def deRHS(s, t):

    n0 = A / (1 + s[M-1] ** 7)
    n = np.hstack([n0, s[M:M + N]])
    N0 = (np.sum(n[1:-2]) + n0 / 2 + n[-1] / 2) * dx  # The integral
    Nt = np.hstack([N0, s[:M]])
    FN = (Nt[:M] - Nt[1:M+1]) / dy
    Fn = (n[:N] - n[1:N+1]) / dx

    return np.hstack([FN, Fn])

# Set default plot parameters
plt.rcParams.update({
    'axes.labelsize': 20,
    'axes.linewidth': 1.2,
    'lines.linewidth': 2.0,
    'patch.linewidth': 0.7,
})

# Steady state analysis
N = np.arange(0, 3, 0.01)
F = 1 / (1 + N ** 7)
blist = [0.2, 0.5, 0.8]

# Plot 1: Red blood cell production
plt.figure(1)
plt.plot(N, F, label='f(U)/A', linewidth=2)
for b in blist:
    plt.plot(N, b * N, label=f'β = {b}', linewidth=2)
plt.xlabel('U')
plt.ylabel('Production rate f(U)/A')
plt.axis([0, 2, 0, 1])
plt.text(1.6, 0.28, 'β = 0.2', fontsize=20)
plt.text(1.6, 0.75, 'β = 0.5', fontsize=20)
plt.text(1.04, 0.8, 'β = 0.8', fontsize=20)
plt.legend()

# Plot 2: U_0 vs. XA
plt.figure(2)
plt.plot(N / F, N, linewidth=2)
plt.ylabel('U_0')
plt.xlabel('XA')
plt.axis([0, 4, 0, 1.4])

# Parameters
d = 7  # Time delay
X = 50  # Lifetime of red blood cells
A = 0.1  # Value for f(0)

# # Plot 3: Stability diagram
x = np.arange(0, 16, 0.01)
f = np.pi * x / ((2 + x) * np.sin(2 * np.pi / (2 + x)))
p = 7
N0p = f / (p - f)
N0 = N0p ** (1 / p)
dA = N0 * (1 + N0p) / x
nnp = np.where(N0p > 0)

plt.figure(3)
plt.plot(x[nnp], dA[nnp], linewidth=2)
plt.xlabel('X/d')
plt.ylabel('d A')
plt.axis([0, 14, 0, 2])
plt.text(5, 1.2, 'Unstable', fontsize=20)
plt.text(5, 0.2, 'Stable', fontsize=20)



# Time-dependent simulation

# Grid
N = 100  # Number of grid points for n
M = 10  # Number of grid points for N
dx = X / (N - 1)  # Grid spacing for n
dy = d / (M - 1)  # Grid spacing for N

# Time span
tspan = np.linspace(0, 500, 500)

# Initial data
s0 = np.hstack([np.ones(M), np.ones(N) * dx / X])

# Solve the ODE
sol = odeint(deRHS, s0, tspan)

n0 = A / (1 + sol[:, M - 1] ** 7)
n = sol[:, M:M + N]
N0 = (np.sum(n[:, :-2], axis=1) + n0 / 2 + n[:, N-1] / 2) * dx  # The integral
    
# Plot the results
plt.figure(4)
plt.plot(tspan, N0, linewidth=2)
plt.xlabel('time (days)')
plt.ylabel('N(t)')



