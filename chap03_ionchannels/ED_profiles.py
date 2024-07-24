
#   -------------------------------------------------------------------
# 
#    Plot profiles from the PNP equations.
# 
#    For Chapter 3, Section 3.3.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

# Set plot parameters
plt.rcParams.update({
    'axes.titlesize': 20,
    'axes.linewidth': 2.0,
    'lines.linewidth': 2.0
})

plt.rcParams['text.usetex'] = True

# Constants
y = np.linspace(0, 1, 100)
ui = 50 / 550
ue = 500 / 550
v = 1
J1 = v * (ui - ue * np.exp(-v)) / (1 - np.exp(-v))
K1 = ui - J1 / v

# Profiles
cprofile = J1 / v + K1 * np.exp(v * y)
vprofile = v - v * y
cprofilelong = ui + (ue - ui) * y
vprofilelong = -(v / np.log(ue / ui)) * np.log(ui / ue + (1 - ui / ue) * y)

# Plot cprofile and vprofile
plt.figure(1)
plt.plot(y, cprofile, label='$u_1$')
plt.plot(y, vprofile, label='$\psi$')
plt.xlabel('y')
plt.legend()

# Plot cprofilelong and vprofilelong
plt.figure(2)
plt.plot(y, cprofilelong, label='$u_1$')
plt.plot(y, vprofilelong, label='$\psi$')
plt.xlabel('y')
plt.legend()

