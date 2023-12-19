#   -------------------------------------------------------------------
# 
#   Plot the steady-state gating variables of the Hodgkin-Huxley model
# 
#    For Chapter 5, Fig. 5.5 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# Set default plot parameters
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.linewidth'] = 2.0
plt.rcParams['lines.linewidth'] = 2.0


Veq = 0  # Original HH formulation
# Veq = -65  # Physiological formulation
Tfact = 1  # Correction factor for temperatures other than 6.3°C
# Tfact = 0.977  # Corresponds to 0°C
# Tfact = 1.085  # Corresponds to 30°C

gnabar = 120.0
gkbar = 36.0
gl = 0.3

Vna = Tfact * (115 + Veq) - Veq
Vk = Tfact * (-12 + Veq) - Veq
Vl = Tfact * (10.5988 + Veq) - Veq

vlist = np.linspace(Vk + 0.01, Vna - 0.01, 100)


def gate_de(V):
    # Calculate the gating functions for V
    AM = 0.1 * (25. - V) / (np.exp(2.5 - 0.1 * V) - 1.)
    BM = 4. * np.exp(-V / 18.)

    AH = 0.07 * np.exp(-V / 20.)
    BH = 1. / (np.exp(0.1 * (30. - V)) + 1.)
    AN = 0.01 * (10. - V) / (np.exp(0.1 * (10. - V)) - 1.)
    BN = 0.125 * np.exp(-V / 80.)

    return np.array([AM, BM, AH, BH, AN, BN])


gates = gate_de(vlist)

taus = np.zeros((3, len(vlist)))
infs = np.zeros((3, len(vlist)))

for j in range(3):
    taus[j, :] = 1. / (gates[2 * j] + gates[2 * j + 1])
    infs[j, :] = gates[2 * j] * taus[j, :]

# Plotting
plt.figure()
plt.plot(vlist, infs[0, :], 'r', label='m_inf(v)')
plt.plot(vlist, infs[1, :], 'b', label='h_inf(v)')
plt.plot(vlist, infs[2, :], 'k', label='n_inf(v)')
plt.legend()
plt.xlabel('Potential (mV)')
plt.show()

plt.figure()
plt.plot(vlist, taus[0, :], 'r', label='tau_m(v)')
plt.plot(vlist, taus[1, :], 'b', label='tau_h(v)')
plt.plot(vlist, taus[2, :], 'g', label='tau_n(v)')
plt.legend()
plt.xlabel('Potential (mV)')
plt.ylabel('ms')
plt.show()
