
#   -------------------------------------------------------------------
# 
#    Compute the steady-state responses and step responses of the simplified
#    modal model.
# 
#    For Chapter 7, Section 7.4.5 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

tau_max = 1000
Ktau = 0.1
tauP = 1
Kc = 0.2
Kh = 0.08
Kp = 0.2
cstim = 0

def get_stuff(c, h, p):
    phi_c = c**4 / (c**4 + Kc**4)
    phi_p = p**2 / (Kp**2 + p**2)
    phi_p_down = Kp**2 / (Kp**2 + p**2)
    h_inf = Kh**4 / (Kh**4 + c**4)
    tauh = tau_max * Ktau**4 / (Ktau**4 + c**4)
    
    beta = phi_p * phi_c * h
    alpha = phi_p_down * (1 - phi_c * h_inf)
    Po = beta / (beta + 0.4 * (beta + alpha))
    
    return tauh, h_inf, Po

def get_c(t):
    global cstim
    return cstim * np.heaviside(t - 2, 1)

def ode_rhs(t, h, p):
    c = get_c(t)
    tauh, h_inf, _ = get_stuff(c, h, p)
    return (h_inf - h) / tauh

# Initialize keep for external plotting
keep = []

# Calculate and plot the steady-state open probability
cc = np.linspace(0.01, 1, 2000)
p_list = [0.2, 0.4, 0.6]
keep.append(cc)

for p in p_list:
    h_inf = Kh**4 / (Kh**4 + cc**4)
    _, _, Po = get_stuff(cc, h_inf, p)
    plt.figure(1)
    plt.plot(cc, Po)
    keep.append(Po)
    
plt.xlim([0, 1])
plt.xlabel('c')
plt.ylim([0, 0.3])
plt.ylabel('P_o')
plt.show()

# Calculate responses to step increases in c
init = [1]
p = 0.4
tspan = np.linspace(0, 30, 2000)
keep.append(tspan)

cstim_list = [0.4, 0.6, 0.8]

for cstim in cstim_list:
    sol = solve_ivp(lambda t, x: ode_rhs(t, x, p), [0, 30], init, t_eval=tspan)
    c = get_c(tspan)
    _, _, Po = get_stuff(c, sol.y[0], p)
    
    plt.figure(2)
    plt.plot(sol.t, Po)
    keep.append(Po)

plt.xlabel('t (s)')
plt.ylabel('P_o')


