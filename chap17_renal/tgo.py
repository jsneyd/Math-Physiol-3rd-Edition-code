
#     -------------------------------------------------------------------
# 
#      Model of tubuloglomerular oscillations.
# 
#      For Chapter 17, Section 17.2 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt


# Define phi function
def phi(x):
    return 1 + K1 * np.tanh(K2 * (Cop - x))


# Define parameters
numx = 100
delx = 1 / numx
delt = 0.001
numtimesteps = 20000
delay = 200
keepend = np.zeros(numtimesteps + delay)
time = np.arange(1, numtimesteps + delay + 1) * delt
c = np.zeros(numx)
c[0] = 1
k = 1
Cop = np.exp(-k)
K1 = 1.5

#  Now solve the PDE by hand (i.e., without using any packaged routines).
#  Take differences in space and time, advance the solution in time, and
#  keep track of the solution at the end of the domain, as it is used for
#  the delay condition inside phi, as well as for the output.

#  The delay is handled in the vector keepend, which we fill up forwards in
#  time (i.e., at timestep i + delay), while using the value at timestep i.
        
gammalist = [3, 2]
for jj, gamma in enumerate(gammalist, start=1):
    K2 = gamma / (K1 * k * np.exp(-k))
    for i in range(numtimesteps):
        keepend[i + delay] = c[numx - 1]
        c[1:numx] = c[1:numx] - delt * k * c[1:numx] - phi(keepend[i]) * (delt / delx) * (c[1:numx] - c[:numx - 1])

    # Plot results
    plt.figure(jj)
    plt.plot(time, keepend)
    plt.xlabel('time')
    plt.ylabel('macula densa [Cl^-]')
    plt.axis([0, 20, 0, 0.7])
    plt.title(f'$\gamma = {gamma}$')

# Stability curves
tbar = np.arange(0, 0.3, 0.001)
for n in range(1, 5):
    w = n * np.pi / (tbar + 1 / 2)
    g = (-1) ** (n + 1) * w / (2 * np.sin(w / 2))
    ndx = np.where(g > 0)
    plt.figure(3)
    plt.plot(tbar[ndx], g[ndx])
    plt.xlabel('$\\tau$')
    plt.ylabel('$\gamma$')
    plt.text(0.1, 19, 'n=4', fontsize=14)
    plt.text(0.17, 15, 'n=3', fontsize=14)
    plt.text(0.2, 5.5, 'n=2', fontsize=14)
    plt.text(0.25, 3, 'n=1', fontsize=14)
    plt.text(0.05, 17, 'unstable', fontsize=14)
    plt.text(0.03, 3, 'stable', fontsize=14)
    plt.axis([0, 0.3, 0, 20])

