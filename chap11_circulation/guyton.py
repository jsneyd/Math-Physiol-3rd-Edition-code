#   -------------------------------------------------------------------
# 
#    Code to solve the autoregulation model of Huntsman/Peskin.
# 
#    For Chapter 11, Section 11.6.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# The data for the three graphs

flow_Pa = np.array([
    [0, 10.256],
    [0, 20.512],
    [0.175, 30.768],
    [0.425, 41.024],
    [0.65, 51.28],
    [0.7875, 61.536],
    [0.875, 71.792],
    [0.9125, 82.048],
    [0.95, 92.30399],
    [0.9875, 102.56],
    [1.025, 112.816],
    [1.05, 123.072],
    [1.075, 133.328],
    [1.1125, 143.584],
    [1.15, 153.84],
    [1.1875, 164.096],
    [1.225, 174.352],
    [1.2875, 184.608],
    [1.35, 194.864],
    [1.45, 205.12],
    [1.6, 215.376],
    [1.8125, 225.632],
    [2.0875, 235.888],
    [2.35, 246.144],
    [2.575, 256.4]
])

flow_metab = np.array([
    [0.829816, 0],
    [0.937826, 0.752941],
    [1.0664, 1.50588],
    [1.23481, 2.25882],
    [1.4571, 3.01176],
    [1.74212, 3.76471],
    [2.0935, 4.51765],
    [2.50963, 5.27059],
    [2.98377, 6.02353],
    [3.50385, 6.77647],
    [4.05269, 7.52941]
])

flow_oxy = np.array([
    [1.02721, 0],
    [1.066, 9.14634],
    [1.11736, 18.2927],
    [1.20706, 27.439],
    [1.35686, 36.5854],
    [1.58445, 45.7317],
    [1.90352, 54.8781],
    [2.3237, 64.0244],
    [2.85057, 73.1707]
])

# Model parameters
O2astar = 104
O2vstar = 40
Pastar = 100
Qstar = 5.6
Mstar = Qstar * (O2astar - O2vstar)
A = 0.1
R0 = 4

# Plot flow against pressure
Pa_model = np.linspace(0, 250, 100)
Qa_model = (1 / (1 + A * O2astar)) * (Mstar * A + Pa_model / R0)/Qstar

plt.figure(1)
plt.plot(flow_Pa[:, 1], flow_Pa[:, 0], label='Data')
plt.plot(Pa_model, Qa_model, label='Model')
plt.xlabel('Arterial pressure (mm Hg)')
plt.ylabel('Blood flow (x normal)')
plt.legend()
plt.show()

# Plot flow against metabolism
M_model = np.linspace(0, 8, 100)
Qa_model2 = (1 / (1 + A * O2astar)) * (M_model * Mstar * A + Pastar / R0)/Qstar

plt.figure(2)
plt.plot(flow_metab[:, 1], flow_metab[:, 0], label='Data')
plt.plot(M_model, Qa_model2, label='Model')
plt.xlabel('Metabolism (x normal)')
plt.ylabel('Blood flow (x normal)')
plt.legend()
plt.show()

# Plot flow against oxygen deficiency
ox_def = np.linspace(0, 50, 100)
O2a = O2astar * (1 - ox_def / 100)
Qa_model3 = (1/(1+A*O2a))*(Mstar*A + Pastar/R0)/Qstar

plt.figure(3)
plt.plot(flow_oxy[:, 1], flow_oxy[:, 0], label='Data')
plt.plot(ox_def, Qa_model3, label='Model')
plt.xlabel('% Arterial oxygen deficiency')
plt.ylabel('Blood flow (x normal)')
plt.legend()
plt.show()
