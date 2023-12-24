#   -------------------------------------------------------------------
# 
#    Use the method of lines to compute the traveling wave in the
#    FitzHugh-Nagumo equations.
# 
#    For Chapter 6, Figure 6.5 and Exercise 6.7 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

# Function to calculate the right-hand side of the ODEs for MOL simulation
def bs_mol(x,t):
    global n, lam

    v = x[:n]
    w = x[n:]
    
    # No-flux boundary condition
    out = np.zeros(2 * n)
    out[0] = lam * (-2 * v[0] + 2 * v[1]) + react(v[0], w[0])
    out[1:n-1] = lam * (v[2:n] - 2 * v[1:n-1] + v[0:n-2]) + react(v[1:n-1], w[1:n-1])
    out[n] = lam * (-2 * v[n-1] + 2 * v[n-2]) + react(v[n-1], w[n-1])

    out[n:] = eps * (v - gamma*w)  # The w equations

    return out

# Function for the reaction term
def react(v, w):
    return v * (alpha - v) * (v - 1) - w + Iapp


# Global parameters
alpha = 0.2
n = 200
L = 150
Iapp = 0
gamma = 0.5
eps = 0.002
diff = 1
delx = L / (n - 1)
lam = diff / (delx ** 2)
x = np.linspace(0, L, n)

# Initial condition
v0 = np.zeros(n)
w0 = np.zeros(n)
v0[:15] = 1

# Time points for output
tspan = np.linspace(0, 200, 10)

# ODE integration
sol = odeint(bs_mol, np.concatenate((v0, w0)), tspan)

# Plot the traveling wave at specific time points
plt.figure(1)
plt.plot(x, sol[1, :n], label=f't = {tspan[1]:.0f}')
plt.plot(x, sol[4, :n], label=f't = {tspan[4]:.0f}')
plt.plot(x, sol[7, :n], label=f't = {tspan[7]:.0f}')
plt.plot(x, sol[9, :n], label=f't = {tspan[9]:.0f}')
plt.xlabel('x')
plt.ylabel('v')
plt.legend()

# Calculate the nullclines
u = np.arange(-0.3, 1.01, 0.01)
w1 = u * (u - alpha) * (1 - u) + Iapp
w2 = u /gamma

# Plot nullclines and add the solution in the phase plane
plt.figure(2)
plt.plot(u, w1, '--', label=r'$\frac{dv}{dt}=0$', linewidth=2)
plt.plot(u, w2, '--', label=r'$\frac{dw}{dt}=0$', linewidth=2)
plt.plot(sol[9, :n],sol[9, n:],linewidth=2)
plt.axis([-0.3, 1, -0.025, 0.15])
plt.xlabel('v')
plt.ylabel('w')
plt.legend()

# plot the 3-d surface
X = np.arange(0, (n - 1) * delx + delx, delx)
fig = plt.figure(3)
ax = fig.add_subplot(111, projection='3d')
X, T = np.meshgrid(X, tspan)
ax.plot_surface(X, T, sol[:, :n], cmap='viridis')
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('v')




