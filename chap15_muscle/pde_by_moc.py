
#   -------------------------------------------------------------------
# 
#    Solve the pde u_y +f(u,t,x) u_x = g(u,t,x) by the method of
#    characteristics.
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

def deRHS(t, s):
    N = len(s) // 2
    x = s[:N]
    u = s[N:]
    
    # Specify the MoC equations
    Fx = u
    Fu = np.zeros(N)
    
    s_prime = np.concatenate((Fx, Fu))
    return s_prime



# Set up the spatial grid
dx = 0.01
x0 = np.arange(0, 1 + dx, dx)
N = len(x0)

# Initial data
u = x0 * (1 - x0)
# Plot initial data
plt.figure()
plt.plot(x0, u, '--', linewidth=2)
plt.xlabel('x', fontsize=20)
plt.ylabel('u', fontsize=20)

# Time parameters
t_end = 3
t_step = 0.25
t_span = np.arange(0, t_end + t_step, t_step)

s0 = np.concatenate((x0, u))
# Solve using solve_ivp
sol = solve_ivp(deRHS, [0, t_end], s0, t_eval=t_span, method='RK45')

# Plot the solution at each time step
for i in range(1, len(sol.t)):
    plt.plot(sol.y[:N, i], sol.y[N:, i], linewidth=2)

plt.show()


