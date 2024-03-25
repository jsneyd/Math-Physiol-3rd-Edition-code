
#   -------------------------------------------------------------------
# 
#   Plot the ventilation-perfusion ratio
# 
#    For Chapter 14, Section 14.5 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

def f(PO):
    Ko2 = 30
    return PO**4 / (Ko2**4 + PO**4)


Pvco2 = 45
Kc = 12
RT = 1.7e4
sigmaco2 = 3.3e-5  # mM/mm Hg
sigmao2 = 1.4e-6
THb = 2.2e-3  # M

Paco2 = np.arange(0.1, 1.01, 0.01) * Pvco2
VbyQ = sigmaco2 * RT * (1 + Kc) * (Pvco2 - Paco2) / Paco2

PiO2 = 150
PvO2 = 40
PaO2 = np.arange(PvO2, 150.1, 0.1)
Ta = sigmao2 * PaO2 + 4 * THb * f(PaO2)
Tv = sigmao2 * PvO2 + 4 * THb * f(PvO2)
VbyQO = RT * (Ta - Tv) / (PiO2 - PaO2)

plt.figure(figsize=(10, 6))
plt.plot(VbyQ, Paco2, label='CO2', color='blue')
plt.plot(VbyQO, PaO2, label='O2', color='red')
plt.xlabel('Ventilation-perfusion ratio')
plt.ylabel('Alveolar partial pressure (mm Hg)')
plt.legend()
plt.grid(True)
plt.xlim(0, 2)
plt.ylim(30, 130)
plt.show()


