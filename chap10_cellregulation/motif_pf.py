
#   -------------------------------------------------------------------
# 
#    Calculate the steady state curves for the Positive feedback motif
# 
#   For Chapter 10, Section 10.1.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

# Set default plot parameters
plt.rcParams.update({
    'axes.labelsize': 20,
    'axes.linewidth': 2.0,
    'lines.linewidth': 2.0,
    'font.size': 20
})

# Parameters
m = 8
n = 4
k1vals = [0.5, 0.2, 1.35]
xlim = [2.5, 4, 2]
Slim = [0.6, 0.3, 2]
k2 = 1

# Plot each curve
for j in range(3):
    k1 = k1vals[j]
    x = np.arange(0, xlim[j] + 0.01, 0.01)

    # Calculate functions
    f = x**m / (1 + x**m)
    y = f / k2
    g = y**n / (1 + y**n)
    S = k1 * x - g

    # Plot
    plt.figure()
    plt.plot(S, x)
    plt.title(f'$k_1 = {k1}$', fontsize=18)
    plt.xlabel('S')
    plt.ylabel('x')
    plt.axis([0, Slim[j], 0, xlim[j]])

