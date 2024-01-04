#   -------------------------------------------------------------------
# 
#    Code for an integrated synapse model that connects the presynaptic
#    voltage to the postsynaptic voltage via calcium and neurotransmitter
#    release.
# 
#    For Chapter 8, Section 8.1.7 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Function to represent the system of ODEs. There are 12 variables
def rhs(sol, t):
    v1, w, o, c, o1, o2, o3, o4, x, y, a, v2 = sol
    R = o1 * o2 * o3 * o4

    k1 = k10 * np.exp(z1 * v1 * FbyRT)
    fac = 2 * v1 * FbyRT
    g = 1 + fac / 2 + fac**2 / 12 - fac**4 / 720 if abs(fac) < 1e-6 else fac / (1 - np.exp(-fac))
    jj = PCa * 2 * F * g * (c - ce * np.exp(-fac))
    ICa = jj * s0 * (o**n)

    tau1 = 1 / (kp1 * c + km1)
    tau2 = 1 / (kp2 * c + km2)
    tau3 = 1 / (kp3 * c + km3)
    tau4 = 1 / (kp4 * c + km4)

    alpha = bigB * np.exp(bigA * v2)
    beta = smallb * np.exp(smalla * v2)

    v1p = 100 * (0.0001 * (v1 - Vr) * (70 - (v1 - Vr)) * ((v1 - Vr) - 7) - w) + Istim / np.cosh(5 * (t % Tperiod) - 1)
    wp = 0.25 * (v1 - Vr - 5 * w)
    op = k1 * (1 - o) - k20 * o
    cp = -ICa * factor - kout * c
    o1p = kp1 * c - o1 / tau1
    o2p = kp2 * c - o2 / tau2
    o3p = kp3 * c - o3 / tau3
    o4p = kp4 * c - o4 / tau4
    xp = -alpha * x + beta * y
    yp = alpha * x + k1m * a * (bigN - x - y) - (beta + k2m) * y
    ap = rho * o1 * o2 * o3 * o4 - ke * a - k1m * a * (bigN - x - y) + k2m * y
    v2p = (1 / Cm) * (-gr * (v2 - Vr) - gs0 * x * (v2 - Vs))

    return [v1p, wp, op, cp, o1p, o2p, o3p, o4p, xp, yp, ap, v2p]

tmax = 15  # this is the length of the simulation
# FHN parameters
Vr=-70             # resting potential 
Tperiod =tmax      # period of stimulus
Istim=100          # stimulus amplitude

# Llinas parameters
k10=2
k20=1
z1=1
n=5

ce=4e4             # units micro molar = 40 mM
PCa=2.0e-2         # mu m ms^{-1}
s0=20

# Voltage-type parameters
F=96490            # Faraday's constant
R=8.315            # Gas constant
Tp=300             # Temperature in degrees Kelvin
FbyRT = F/(R*Tp*1000) # units of mV^(-1)

# Calcium parameters. Radius in microns. Factor converts current to calcium flux.
factor=(3/20)*(1/(2*F))
kout = 3           # calcium withdrawal rate

# Calcium-stimulated secretion parameters for o_j dynamics
kp1=3.75e-3
km1=4e-4     
kp2=2.5e-3 
km2=1e-3   
kp3=5e-4
km3=0.1
kp4=7.5e-3
km4=10
rho = 3e6

# Magleby parameters
bigA=0.008
bigB=1.43
bigN=10            #dimensionless
smalla=0.00315 
smallb=1
k1m=1000           # ms^(-1) 
k2m=500            # mM ms^(-1)
ke=10              # rate of ACh degradation

# Postsynaptic voltage parameters
Vs=-15 
gr=10
gs0 = 10
Cm=1



# Initial conditions
init = [-70, 0, 0.11, 1.4e-6, 0.14, 0.05, 0.0015, 0, 0, 0, 0, -70]

# Time points
t = np.arange(0, 15, 0.01)

# Solve ODEs
sol = odeint(rhs, init, t)

# Plot results
plt.figure(figsize=(10, 8))
plt.subplot(3, 3, 1)
plt.plot(t, sol[:, 0])
plt.xlabel('t (ms)')
plt.ylabel('V1 (ms)')

plt.subplot(3, 3, 2)
plt.plot(t, sol[:, 2])
plt.ylabel('o')
plt.xlabel('t (ms)')

plt.subplot(3, 3, 3)
plt.plot(t, sol[:, 3])
plt.ylabel('c (uM)')
plt.xlabel('t (ms)')

plt.subplot(3, 3, 4)
plt.plot(t, sol[:, 4], label='o1')
plt.plot(t, sol[:, 5], label='o2')
plt.plot(t, sol[:, 6], label='o3')
plt.plot(t, sol[:, 7], label='o4')
plt.legend()
plt.xlabel('t (ms)')

plt.subplot(3, 3, 5)
plt.plot(t, sol[:, 4] * sol[:, 5] * sol[:, 6] * sol[:, 7])
plt.ylabel('PR')
plt.xlabel('t (ms)')

plt.subplot(3, 3, 6)
plt.plot(t, sol[:, 10])
plt.xlabel('t (ms)')
plt.ylabel('a (mM)')

plt.subplot(3, 3, 7)
plt.plot(t, sol[:, 8], label='x')
plt.plot(t, sol[:, 9], label='y')
plt.legend()
plt.xlabel('t (ms)')

plt.subplot(3, 3, 8)
plt.plot(t, sol[:, 0], label='presynaptic')
plt.plot(t, sol[:, 11], label='postsynaptic')
plt.legend()
plt.xlabel('t (ms)')
plt.ylabel('Membrane Potential (ms)')

plt.tight_layout()
plt.show()
