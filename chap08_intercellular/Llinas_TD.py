#   -------------------------------------------------------------------
# 
#    Program to solve the Llinas model of synaptic suppression, with a 
#    time-dependent voltage input. The voltage is computed as a solution of
#    the FHN equations.
# 
#    For Chapter 8, Section 8.1.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Define the system of ODEs
def rhs(s, t):
    v1, w, op = s
    
    # FHN odes for the presynaptic voltage
    v1p = 100 * (0.0001 * (v1 - Vr) * (70 - (v1 - Vr)) * ((v1 - Vr) - 7) - w)
    wp = 0.25 * (v1 - Vr - 5 * w)
    
    # Linas ode for open probability 
    k1 = k10 * np.exp(FbyRT * z1 * v1)
    op = k1 * (1 - op) - k2 * op
    
    return [v1p, wp, op]

def fCa(V,oh):
    phi=2*FbyRT*V;
    jj=PCa*2*F*phi*(ci-ce*np.exp(-phi))/(1-np.exp(-phi));
    ICa =(s0/n)*jj*(oh**n);
    return ICa

# Constants
k10 = 2
k2 = 1
z1 = 1
s0 = 100
PCa = 0.00001
F = 96490
R = 8.315
Tp = 300
ci = 0
ce = 40
n = 5
Vss = -70
FbyRT = F / (R * Tp * 1000)
k1 = k10 * np.exp(FbyRT * z1 * Vss)
ohss = k1 / (k1 + k2)
Vr = -70

# Time vector
t = np.arange(0, 3.01, 0.01)

# Initial conditions
init = [-50, 0, 0.11]


# Solve the system of ODEs
sol = odeint(rhs, init, t)
V = sol[:, 0]
oh = sol[:, 2]

# Compute calcium current
current = fCa(V, oh)

# Plot the solutions
plt.figure(1)
plt.plot(t, V, label='V', linewidth=2)
plt.xlabel('time (ms)')
plt.ylabel('V (mV)')

plt.figure(2)
fig, ax1 = plt.subplots()
ax2 = ax1.twinx() 
ax1.plot(t, oh,  'b', label='open probability')
ax1.tick_params(axis ='y', colors = 'blue')
ax1.set_xlabel('time (ms)')
ax1.set_ylabel('open probability',color='blue')
ax1.spines['left'].set_color('blue')

ax2.plot(t, current, 'r', label=r'$I_{\rm Ca}$', linewidth=2)
ax2.set_ylabel(r'$I_{\rm Ca}$ (pA/($\mu$m)$^2$)',color='red')
ax2.tick_params(axis ='y', colors = 'red')
ax2.spines['right'].set_color('red')


