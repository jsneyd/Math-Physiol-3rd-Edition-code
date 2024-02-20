
#   -------------------------------------------------------------------
# 
#    Simulation of the 6-variable generic model of the cell cycle
# 
#    For Chapter 10, Section 10.4.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the differential equations
def ccgen5odes(y, t):
    cycB, cdh, cdcT, cdcA, IEP, m = y
    
    dy1 = k1 * m - (k2p + k2pp * cdh) * cycB
    dy2 = (k3p + k3pp * cdcA) * (1 - cdh) / (J3 + 1 - cdh) - k4 * cycB * cdh / (J4 + cdh)
    dy3 = k5p + k5pp * (cycB ** 4) / (J5 ** 4 + (cycB ** 4)) - k6 * cdcT
    dy4 = k7 * IEP * (cdcT - cdcA) / (J7 + cdcT - cdcA) - k8 * cdcA / (J8 + cdcA) - k6 * cdcA
    dy5 = k9 * cycB * (1 - IEP) - k10 * IEP
    dy6 = mu * m
    
    return [dy1, dy2, dy3, dy4, dy5, dy6]

# Parameters
k1 = 0.04
k2p = 0.04
k2pp = 1.0
k2ppp = 1.0
k3p = 1.0
k3pp = 10.0
k4p = 2.0
k4 = 35.0
k5p = 0.005
k5pp = 0.2
k6 = 0.1
k7 = 1.0
k8 = 0.5
J3 = 0.04
J4 = 0.04
J5 = 0.3
J7 = 1.0e-3
J8 = 1.0e-3
k9 = 0.1
k10 = 0.02
mu = 0.005

npts=50
tend=1
tfinal=300

#  variables are (in this order)
#  cycB;
#  cdh ;
#  cdcT ;
#  cdcA ;
#  IEP ;
#  m ;

# Initial conditions
y0 = [0.0204, 0.9148, 0.05, 0.0001, 0.062, 0.4943]

# Time points
t = np.linspace(0, 300, 500)

# Solve the differential equations
keep = np.zeros((tfinal, 7))  # Initialize array to store results
keep[:, 0] = np.arange(tfinal) * tend  # Time points

y0 = [0.0204, 0.9148, 0.05, 0.0001, 0.062, 0.4943]  # Initial conditions

test = 0  # Initialize test variable
for loop in range(tfinal):
    tspan = np.linspace(0, tend, npts)
    sol = odeint(ccgen5odes, y0, tspan)
    # Store the solutions
    keep[loop, 1:7] = sol[-1, :]

    # Check condition for CycB concentration
    if sol[-1, 0] > 0.4:
        test = 1
    # Reset mass for cell division if condition is met
    if sol[-1, 0] < 0.05 and test == 1:
        sol[-1, 5] /= 2.0
        test = 0

    # Update initial conditions for the next iteration
    y0 = sol[-1, :]


# Plot the results
plt.figure()
plt.plot(keep[:,0],keep[:,1], label='CycB')
plt.plot(keep[:,0],keep[:,2], label='Cdh')
plt.plot(keep[:,0],keep[:,3], label='Cdc_T')
plt.plot(keep[:,0],keep[:,6], label='m')
plt.title('6-Variable Generic Model of the Cell Cycle')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

