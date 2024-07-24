
#   -------------------------------------------------------------------
# 
#    Comparing the GHK and linear IV curves.
# 
#    For Chapter 3, Section 3.1.1 of
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
FRT = 1 / 25.8
nai = 50
nae = 437
vna = (1 / FRT) * np.log(nae / nai)
gna = 0.01
pnaF = -gna * vna / (nai - nae)  # to guarantee the curves meet at V=0

ki = 397
ke = 20
gk = 0.367
vk = (1 / FRT) * np.log(ke / ki)
pkF = -gk * vk / (ki - ke)  # this is fixed at value determined by normal kae

# Voltage range
x = np.linspace(-40, 100, 100)
# Linear IV curve for Na+
linIV = gna * (x - vna)
# GHK IV curve for Na+
GHKIV = pnaF * FRT * x * (nai - nae * np.exp(-x * FRT)) / (1 - np.exp(-x * FRT))

# Potassium concentration range
x1 = np.linspace(10, 50, 100)
# Linear Vr curve
linVr = (gna * vna + gk * (1 / FRT) * np.log(x1 / ki)) / (gna + gk)
# GHK Vr curve
GHKVr = (1 / FRT) * np.log((pnaF * nae + pkF * x1) / (pnaF * nai + pkF * ki))

# Plot IV curves
plt.figure(1)
plt.plot(x, linIV, label='linear')
plt.plot(x, GHKIV, label='GHK')
plt.legend()
plt.xlabel('$V$ (mV)')
plt.ylabel('$I_{Na}$')
plt.title('IV Curves')
plt.show()

# Plot Vr curves
plt.figure(2)
plt.plot(x1, linVr, label='linear')
plt.plot(x1, GHKVr, label='GHK')
plt.legend()
plt.xlabel('[K$^+$]$_e$ (mM)')
plt.ylabel('$V_r$ (mV)')
plt.title('$V_r$ Curves')
plt.show()
