

#   -------------------------------------------------------------------
# 
#    To solve u_t + J_x= 0, with MoL and upwinding (written in conservation
#    form, J = v*u).
# 
#    For Chapter 15, Section 15.12.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def deRHS(t, u):
    
    # Evaluate v at x_{j-1/2}
    #xjmh = np.concatenate(([x - dx / 2], [x[N-1] + dx / 2]))
    #ujmh = (np.concatenate(([u], [0])) + np.concatenate(([0], [u]))) / 2

    # Evaluate v_{j-1/2} = v(xjmh, t, ujmh)
    vjmh = -2 * np.sin(3 * t)  # v as a function of time
    
    # Upwinding
    Jmh = vjmh * ((vjmh > 0) * np.concatenate(([0], [u][0])) 
               + (vjmh < 0) * np.concatenate(([u][0], [0])))
    Fu = (Jmh[:-1] - Jmh[1:]) / dx  # Finite difference in x

    rxn = ((x > 0) & (x < 1)) * (1 - u) - 0.5 * u
    s_prime = Fu + rxn
    
    return s_prime


# Spatial grid
dx = 0.01
x = np.arange(-3, 3 + dx, dx)  # grid points
N = len(x)  # number of points in grid

# Initial data
u = np.exp(-x**2)

# Time parameters
t_end = 10
t_step = 0.1
t_span = np.linspace(0,t_end,101)

# Solve the MoL differential equations
u0 = u  # initial condition as a 1D array

# Solve using solve_ivp
sol = solve_ivp(deRHS, [0, t_end], u0, t_eval=t_span, method='RK45')

# Plotting results
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#T, X = np.meshgrid(sol.t, x)
for i in range(len(sol.t)):
    ax.plot(x, sol.y[:, i], sol.t[i])


ax.set_xlabel('x')
ax.set_ylabel('u')
ax.set_zlabel('t')
ax.view_init(vertical_axis='y')
plt.show()



