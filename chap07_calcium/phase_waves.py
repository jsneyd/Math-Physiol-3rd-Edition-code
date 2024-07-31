
#   -------------------------------------------------------------------
# 
#    Use the method of lines to compute a traveling calcium wave. This
#    version includes IP3 diffusion, and is used to demonstrate phase waves.
# 
#    For Chapter 7, Section 7.9 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Function to compute the right-hand side of the PDE
def pde_rhs(t, s):
    P = s[:N]
    C = s[N:2*N]
    
    scu = du / h**2
    scv = dv / h**2
    
    out = coscrhs(t, s)
    
    Fp = out[:N]
    Fc = out[N:2*N]
    Fh = out[2*N:3*N]
    
    FP = scu * (-sc * P + np.concatenate(([0], P[:-1])) + np.concatenate((P[1:], [0]))) + Fp
    FC = scv * (-sc * C + np.concatenate(([0], C[:-1])) + np.concatenate((C[1:], [0]))) + Fc
    
    s_prime = np.concatenate([FP, FC, Fh])
    
    return s_prime

# Function to compute the right-hand side of the ODE
def coscrhs(t, s):
    P = s[:N]
    C = s[N:2*N]
    h = s[2*N:3*N]
    ce = gamma * (ct - C)
    
    phi_c = C**4 / (C**4 + Kc**4)
    phi_p = P**2 / (Kp**2 + P**2)
    phi_p_down = Kp**2 / (Kp**2 + P**2)
    h_inf = Kh**4 / (Kh**4 + C**4)
    tauh = tau_max * Ktau**4 / (Ktau**4 + C**4)
    
    beta = phi_p * phi_c * h
    alpha = phi_p_down * (1 - phi_c * h_inf)
    Po = beta / (beta + 0.4 * (beta + alpha))
    serca = Vserca * (C**2 - Kbar * ce**2) / (C**2 + kserca**2)
    
    Fp = np.zeros(N)
    Fc = kf * Po * (ce - C) - serca
    Fh = (h_inf - h) / tauh
    
    out = np.concatenate([Fp, Fc, Fh])
    
    return out

plt.rcParams.update({
    'axes.titlesize': 20,
    'axes.linewidth': 2.0,
    'lines.linewidth': 2.0
})
plt.rcParams['text.usetex'] = True


# Parameters
ct = 2
tau_max = 1000
Ktau = 0.1
tauP = 1
Kc = 0.2
Kh = 0.08
kf = 10
Vserca = 0.9
kserca = 0.2
Kbar = 0.00001957
Kp = 0.2
gamma = 5.5
N = 300
L = 10

du = 0.3    # IP3 diffusion coefficient
dv = 0.001  # Ca diffusion coefficient


h = L / N
sc = np.concatenate(([1], 2 * np.ones(N - 2), [1]))

X = h * np.arange(1, N + 1)

# Initial conditions
P0 = 1.5 * (X < 1).astype(float)
C0 = 0.1 * np.ones(N)
h0 = 0.5 * np.ones(N)
init = np.concatenate([P0, C0, h0])

# Time integration
tstep = 0.5
t_end = 150
tspan = np.arange(0, t_end + tstep, tstep)

sol = solve_ivp(lambda t, x: pde_rhs(t, x), [0, t_end], init, t_eval=tspan, max_step=1)

P = sol.y[:N]
C = sol.y[N:2*N]
h = sol.y[2*N:3*N]

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(X, P[:, -1], label='IP$_3$')
plt.plot(X, C[:, -1], label='Ca$^{++}$')
plt.legend(fontsize=20)
plt.xlabel('x', fontsize=20)
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(sol.t, C[int(N/2)], label='Ca$^{++}$(L/2)')
plt.xlabel('time', fontsize=20)
plt.ylabel('Ca$^{++}$(L/2)', fontsize=20)
plt.show()

plt.figure(figsize=(10, 6))
plt.contour(X, sol.t, C.T, linewidths=2)
plt.xlabel('x', fontsize=20)
plt.ylabel('t', fontsize=20)
plt.show()
