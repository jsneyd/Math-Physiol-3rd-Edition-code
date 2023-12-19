
#   -------------------------------------------------------------------
#  
#    Python code to calculate the fundamental solution and several Green's
#    functions for the cable equation on a finite domain. The Green's 
#    functions are calculated
#    for sealed-end and closed-circuit ends (i.e., Neumann and Dirichlet BCs)
#    and are calculated by a sum of fundamental solutions and by Fourier
#    series. So there are several different curves.
#  
#    For Chapter 4 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

def fund(t, x, xi, tau):
    d1 = 1. / np.sqrt(4 * np.pi * (t - tau))
    d2 = np.exp(-(x - xi) ** 2 / (4 * (t - tau)))
    d3 = np.exp(-(t - tau))
    return d1 * d2 * d3

def baseF(t, x, xi, tau, n, L):
    d1 = np.cos(n * np.pi * xi / L) * np.cos(n * np.pi * x / L)
    d2 = np.exp(-(1 + (n * np.pi / L) ** 2) * (t - tau))
    return d1 * d2 / L

def baseFSC(t, x, xi, tau, n, L):
    d1 = np.sin(n * np.pi * xi / L) * np.sin(n * np.pi * x / L)
    d2 = np.exp(-(1 + (n * np.pi / L) ** 2) * (t - tau))
    return d1 * d2 / L

def GMISE(T, X, xi, tau, L):
    return (fund(T, X, xi, tau) + fund(T, X, -xi, tau) + fund(T, X, 2 * L - xi, tau) +
            fund(T, X, 2 * L + xi, tau) + fund(T, X, -2 * L + xi, tau))

def GMISC(T, X, xi, tau, L):
    return (fund(T, X, xi, tau) - fund(T, X, -xi, tau) + fund(T, X, 2 * L + xi, tau) -
            fund(T, X, -2 * L - xi, tau) + fund(T, X, 3 * L + xi, tau) - fund(T, X, -3 * L - xi, tau) -
            fund(T, X, 2 * L - xi, tau) + fund(T, X, -2 * L + xi, tau) - fund(T, X, 3 * L - xi, tau) +
            fund(T, X, -3 * L + xi, tau) - fund(T, X, 4 * L - xi, tau))

def GFSSE(t, x, xi, tau, L):
    out = baseF(t, x, xi, tau, 0, L)
    for i in range(1, 21):
        out += baseF(t, x, xi, tau, i, L) + baseF(t, x, xi, tau, -i, L)
    return out

def GFSSC(t, x, xi, tau, L):
    out = baseFSC(t, x, xi, tau, 0, L)
    for i in range(1, 21):
        out += baseFSC(t, x, xi, tau, i, L) + baseFSC(t, x, xi, tau, -i, L)
    return out


# End of the functions. Now get stuff done. I mean, not real work. 
# I'm a professor, I don't do that, but kind of stuff that isn't just
# defining functions.

L = 1
xi = 0.7 * L
tau = 0

# Plot fundamental solution on the infinite domain
T, X = np.meshgrid(np.arange(tau, tau + 0.2, 0.005), np.arange(xi - 0.8, xi + 0.8, 0.01))
Z = fund(T, X, xi, tau)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot_surface(T, X, Z, cmap='viridis')
ax.set_xlabel('T')
ax.set_ylabel('X')
ax.set_zlabel('V')
ax.set_title('Fundamental solution on the infinite domain')

# Plot Green's function on finite domain with sealed ends
T, X = np.meshgrid(np.arange(tau, tau + 0.2, 0.001), np.arange(0, L, L / 200, dtype=float))
Z1 = GMISE(T, X, xi, tau, L)
Z2 = GFSSE(T, X, xi, tau, L)
fig = plt.figure()
ax1 = fig.add_subplot(projection='3d')
ax1.plot_surface(T, X, Z1, cmap='viridis')
ax1.set_xlabel('T')
ax1.set_ylabel('X')
ax1.set_zlabel('V')
ax1.set_title('Method of images (Sealed Ends)')
plt.show()

fig = plt.figure()
ax2 = fig.add_subplot(projection='3d')
ax2.plot_surface(T, X, Z2, cmap='viridis')
ax2.set_xlabel('T')
ax2.set_ylabel('X')
ax2.set_zlabel('V')
ax2.set_title('Fourier expansion (Sealed Ends)')
plt.show()


# Plot for fixed T values
tf = [0.001, 0.01, 0.05, 0.1]
xf = np.linspace(0, L, 500)
Vf1 = GMISE(tf[0], xf, xi, tau, L)
Vf2 = GMISE(tf[1], xf, xi, tau, L)
Vf3 = GMISE(tf[2], xf, xi, tau, L)
Vf1F = GFSSE(tf[0], xf, xi, tau, L)
Vf2F = GFSSE(tf[1], xf, xi, tau, L)
Vf3F = GFSSE(tf[2], xf, xi, tau, L)
plt.figure()
plt.plot(xf, Vf1, 'k', xf, Vf2, 'k', xf, Vf3, 'k', linewidth=2)
plt.plot(xf, Vf1F, 'r--', xf, Vf2F, 'b--', xf, Vf3F, 'g--', linewidth=2)
plt.xlabel('T')
plt.ylabel('V')
plt.legend(['Method of Images', '_', '_', 'Fourier expansion: t=0.001', 't=0.01', 't=0.05'], loc='best')
plt.title('Comparison for fixed T values')
plt.show()


# Plot Green's function on finite domain with short-circuit ends
T, X = np.meshgrid(np.arange(tau, tau + 0.2, 0.005), np.arange(0, L, L / 200, dtype=float))
Z1 = GMISC(T, X, xi, tau, L)
Z2 = GFSSC(T, X, xi, tau, L)
fig = plt.figure()
ax1 = fig.add_subplot(projection='3d')
ax1.plot_surface(T, X, Z1, cmap='viridis')
ax1.set_xlabel('T')
ax1.set_ylabel('X')
ax1.set_zlabel('V')
ax1.set_title('Method of images (Short-Circuit Ends)')

fig = plt.figure()
ax2 = fig.add_subplot(projection='3d')
ax2.plot_surface(T, X, Z2, cmap='viridis')
ax2.set_xlabel('T')
ax2.set_ylabel('X')
ax2.set_zlabel('V')
ax2.set_title('Fourier expansion (Short-Circuit Ends)')
plt.show()

# Plot for fixed T values
tf = [0.001, 0.01, 0.05, 0.1]
xf = np.linspace(0, L, 500)
Vf1 = GMISC(tf[0], xf, xi, tau, L)
Vf2 = GMISC(tf[1], xf, xi, tau, L)
Vf3 = GMISC(tf[2], xf, xi, tau, L)
Vf1F = GFSSC(tf[0], xf, xi, tau, L)
Vf2F = GFSSC(tf[1], xf, xi, tau, L)
Vf3F = GFSSC(tf[2], xf, xi, tau, L)
fig = plt.figure()
plt.plot(xf, Vf1, 'k', xf, Vf2, 'k', xf, Vf3, 'k', linewidth=2)
plt.plot(xf, Vf1F, 'r--', xf, Vf2F, 'b--', xf, Vf3F, 'g--', linewidth=2)
plt.xlabel('X')
plt.ylabel('V')
plt.legend(['Method of images', '_', '_', 'Fourier expansion: t=0.001', 't=0.01', 't=0.05'])
plt.title('Comparison for fixed T values')
plt.show()
