#   -------------------------------------------------------------------
# 
#    Program to solve the Llinas model of synaptic suppression.
# 
#    For Chapter 8, Section 8.1.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# Calculate oh analytically
def prob(t, V):
    k1 = k10 * np.exp(FbyRT * z1 * V)
    oh = k1 / (k1 + k2) * (1 - np.exp(-(k1 + k2) * t)) + ohss * np.exp(-(k1 + k2) * t)
    return oh


# calculate the currents
def fCa(V, oh):
    phi = 2 * FbyRT * V
    jj = PCa * 2 * F * phi * (ci - ce * np.exp(-phi)) / (1 - np.exp(-phi))
    ICa = (s0 / n) * jj * (oh ** n)
    return ICa



# Calculate oh using the differential equation
def probde(t, V):
    k1 = k10 * np.exp(FbyRT * z1 * Vss)
    ohss = k1 / (k1 + k2)
    sinit = ohss
    tspan = t

    S = odeint(rhs, sinit, tspan, args=(V,))
    oh = S[:, 0]
    return oh


# Define the right-hand side of the differential equation
def rhs(s, t, V):
    Vt = np.where(t < tsw, V, Vss)
    k1 = k10 * np.exp(FbyRT * z1 * Vt)
    return k1 * (1 - s) - k2 * s


k10 = 2
k2 = 1
z1 = 1
s0 = 100
PCa = 0.00001
F = 96490
R = 8.315
Tp = 300
ci = 0.00001
ce = 40
n = 5
Vss = -70
FbyRT = F / (R * Tp * 1000)  # units of mV^(-1)
k1 = k10 * np.exp(FbyRT * z1 * Vss)
ohss = k1 / (k1 + k2)
tsw = 2.5

# Time vector
t = np.arange(0, 5.01, 0.01)

# Voltage list
Vlist = [-40, -10, 20, 50, 120]

# SOlution for different voltage step increases
plt.figure(1)
ICa_save = np.zeros((len(Vlist), len(t)))
for j, V in enumerate(Vlist):
    oh = prob(t, V)
    ICa = fCa(V, oh)
    
    plt.plot(t, -ICa, linewidth=2)
    plt.xlabel('time (ms)')
    plt.ylabel('-I_{Ca}(pA/(\mum)^2)')
    plt.legend(['V=-40mV', 'V=-10mV', 'V=20mV', 'V=50mV', 'V=120mV'])



# Steady-state curves
x = np.linspace(-70, 70, 1000)
phi_x = 2 * FbyRT * x
Jss = PCa * 2 * F * phi_x * (ci - ce * np.exp(-phi_x)) / (1 - np.exp(-phi_x))
Oss = k10 * np.exp(FbyRT * z1 * x) * (1 / (k10 * np.exp(FbyRT * z1 * x) + k2))
ICass = Jss * (s0 / n) * (Oss ** n)

# Plot the steady-state curves
plt.figure(2)
fig, ax1 = plt.subplots()
ax2 = ax1.twinx() 
ax1.plot(x, ICass, 'r', x, Jss,'b', linewidth=2)
ax2.plot(x, Oss, 'g', linewidth=2)
ax1.set_xlabel('V (mV)')
ax1.set_ylabel(r'$I_{Ca}$')
ax2.set_ylabel('open probability',color='green')
ax2.tick_params(axis ='y', colors = 'green')
ax2.spines['right'].set_color('green')
plt.legend()
plt.show()


# Plot the responses to a step increase in V followed by a step decrease
for j, V in enumerate(Vlist):
    Vt = np.where(t < tsw, V, Vss)
    oh = probde(t, V)
    ICa = fCa(Vt, oh)

    plt.figure(3)
    plt.semilogy(t, -ICa, linewidth=2)
    plt.xlabel('time (ms)')
    plt.ylabel(r'$-I_{Ca}$ (pA/($\mu$ m)$^2)$')
    plt.legend(['V=-40mV', 'V=-10mV', 'V=20mV', 'V=50mV', 'V=120mV'])

    plt.figure(4)
    plt.plot(t, oh)
    plt.xlabel('time (ms)')
    plt.ylabel('open probability')

plt.figure(3)
plt.ylim([0.1,1e4])
plt.legend(['V=-40mV', 'V=-10mV', 'V=20mV', 'V=50mV', 'V=120mV'])
plt.figure(4)
plt.legend(['V=-40mV', 'V=-10mV', 'V=20mV', 'V=50mV', 'V=120mV'])




