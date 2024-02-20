import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def ddeRHS(t, Y):
    P, H = Y
    P_delay = P(t - tau)
    dPdt = -P / (ca * R) + (Vs / ca) * H
    dHdt = alpha * gs(P_delay) - beta * (1 - gs(P))
    return [dPdt, dHdt]

def gs(x):
    return mu**n / (mu**n + x**n)

# Global Parameters
n = 7
tau = 2
ca = 1.55
R = 1.05
Vs = 67.9 
alpha = 0.84
beta = 1.17
mu = 93.0

# Main Simulation
t_span = (0, 200)
init = [88.71759, 0.8026246]  # Initial conditions


# Solve DDE system
sol = solve_ivp(ddeRHS, t_span, init, t_eval=np.linspace(*t_span, 1000), method='RK45')

# Plot results
plt.figure(figsize=(10, 5))
plt.plot(sol.t, sol.y[0], 'r', label='P')
plt.plot(sol.t, sol.y[1], 'b', label='H')
plt.xlabel('t')
plt.legend()
