
#   -------------------------------------------------------------------
# 
#    Ectopic focus oscillatory waves
# 
#    For Chapter 12, Section 12.4.4 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Define the function for the ODE (Method of Lines)
def deRHS(t, u):
    v = u[:N]
    w = u[N:]

    vp = Dscal * np.matmul(A, v) + f(v) - w
    wp = eps * (v - gam * w - alps)

    return np.concatenate((vp, wp))

# Define the function for f(v)
def f(u):
    return 10 * u * (1 - u) * (u - 0.5)

# Define the parameters
eps = 0.1
gam = 0.1
alf = 0.0
b = 0.5
scale = 1.5
N = 200
L = 20

dx = L / N
x = np.arange(1, N + 1) * dx
alps = alf + (b - alf) * np.exp(-((x / scale) ** 2))

Dl = 1
Dscal = Dl / (dx ** 2)
# impose a no-flux Neumann boundary condition at both ends
A = -2 * np.diag(np.ones(N)) + np.diag([2] + [1]*(N - 2), 1) + np.diag([1]*(N - 2) + [2], -1)

V0 = np.zeros(N)
W0 = np.zeros(N)
tstep = 0.5
t_end = 150
tspan = (0, t_end)
s0 = np.concatenate((V0, W0))

# Solve the ODE
sol = solve_ivp(deRHS, tspan, s0, t_eval=np.arange(0, t_end, tstep),method='Radau')

# Plot the solution
plt.figure()
ax = plt.axes(projection='3d')
X, T = np.meshgrid(x, sol.t)
ax.plot_surface(X, T, np.transpose(sol.y[:N,:]))
plt.xlabel('x')
plt.ylabel('t')

