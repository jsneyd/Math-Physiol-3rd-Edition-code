#   -------------------------------------------------------------------
# 
#    This code sets n and h to be constant and looks at the nullclines in
#    the fast Hodgkin-Huxley model
# 
#    For Chapter 5, Fig. 5.8 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

# Set default plot parameters
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['lines.linewidth'] = 2.0


# Define global parameters
Veq = 0  # Original HH formulation
# p.Veq = -65  # Physiological formulation
Tfact = 1  # Correction factor for temperatures other than 6.3°C
# Tfact = 0.977  # Corresponds to 0°C
# Tfact = 1.085  # Corresponds to 30°C

gnabar = 120.0
gkbar = 36.0
gl = 0.3
Vna0 = 115
Vk0 = -12
Vl0 = 10.5988

Vna = Tfact * (Vna0 + Veq) - Veq
Vk = Tfact * (Vk0 + Veq) - Veq
Vl = Tfact * (Vl0 + Veq) - Veq

Inp = 0  # This is the amplitude of the initial input stimulus
# p.Inp = 3.47  # This is the amplitude of the initial input stimulus

# Iapp is the amplitude of the steady current stimulus
Iapp = 0  # Pick p.Iapp < 0 to make things easier to see
BCL = 39
# p.BCL = 0

# Assign initial data
V = Veq
m = 0.0529
h0set = [0.5961, 0.4, 0.2, 0.1]
n0set = [0.3177, 0.5, 0.7, 0.8]

# first plot the dm/dt=0 nullcline
Vm = np.arange(-20, 120, 0.1)
AM = 0.1 * (25. - Vm) / (np.exp(0.1 * (25. - Vm)) - 1.)
BM = 4. * np.exp(-Vm / 18)
Minf = AM / (AM + BM)
plt.plot(Vm, Minf, 'b', label='Nullcline m_inf(v)')
plt.xlabel('V')
plt.ylabel('m')
plt.axis([-10, 120, 0, 1])
plt.title(f'I_app = {Iapp:.2f}', fontsize=18)

# now plot a selection of dv/dt=0 nullclines
for j in range(len(n0set)):
    n0 = n0set[j]
    h0 = h0set[j]
    M = np.arange(0, 1, 0.001)
    gk = gkbar * n0 ** 4
    gna = gnabar * M ** 3 * h0

    Vs = (gna * Vna + gk * Vk + gl * Vl + Iapp) / (gna + gk + gl)
    plt.plot(Vs, M, 'g', label='Nullcline dv/dt=0')
    plt.legend(['dm/dt=0','dv/dt=0','_','_','_'])

