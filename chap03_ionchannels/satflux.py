#   -------------------------------------------------------------------
# 
#    Flux in a saturating one-ion model.
# 
#    For Chapter 3, Section 3.4.2 of
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
k0 = 1
k1 = 1
km1 = 1
km2 = 1

alpha1 = (km1 + k1) / (k0 * k1)
beta1 = 1 / k1
gamma1 = km2 / (k0 * k1)

ce = 1

# Intracellular concentration
ci = np.linspace(0, 10, 100)

# Flux
J = (ci - ce * (km1 * km2 / (k0 * k1))) / (alpha1 + beta1 * ci + gamma1 * ce)

# Plot
plt.plot(ci, J)
plt.xlabel('$c_i$')
plt.ylabel('$J$')
plt.title('Flux $J$ as a function of intracellular concentration $c_i$')

