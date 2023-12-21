 
#   -------------------------------------------------------------------
# 
#   This is an ode integrator for the FitzHugh-Nagumo equations.
# 
#    For Chapter 5 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def fhnrhs(t, x):
    v, w = x
    dvdt = (v * (alpha - v) * (v - 1) - w + Iapp) / eps
    dwdt = v - gamma * w
    return [dvdt, dwdt]

eps = 0.01
alpha = 0.1
gamma = 0.5
Iapp = 0

t_span = (0, 4)
tout = np.linspace(0,4,500)
initial_conditions = [0.2, 0]

sol = solve_ivp(fhnrhs, t_span, initial_conditions,t_eval=tout)

# Plotting
plt.figure(1)
plt.plot(sol.t, sol.y[0], 'r', label='v')
plt.plot(sol.t, sol.y[1], 'b', label='w')
plt.xlabel('Time')
plt.ylabel('v')
plt.legend()

plt.figure(2)
v = np.linspace(-0.4, 1.4, 100)
w1 = v * (alpha - v) * (v - 1) + Iapp
w2 = v / gamma
plt.plot(v, w1, 'g--', label='dv/dt=0')
plt.plot(v, w2, 'b--', label='dw/dt=0')
plt.plot(sol.y[0], sol.y[1], 'r', label='Trajectory')
plt.ylim([-0.1, 1.2 * max(w1)])
plt.xlabel('v')
plt.ylabel('w')
plt.legend()

plt.show()
