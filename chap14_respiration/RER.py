
#   -------------------------------------------------------------------
# 
#   Plot the respiratory exchange ratio.
# 
#    For Chapter 14, Sections 14.5.1 and 14.5.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

def hemo(PO):
    return PO**4 / (KO2**4 + PO**4)

THb = 2.2e-3  # M
sigmaO2 = 1.4e-6
sigmaCO2 = 3.3e-5
PvO2 = 40.0
PvCO2 = 45.0
PiO2 = 150.0
RT = 1.7e4
Z0 = 2.2e-3
KO2 = 30 * sigmaO2
Kc = 12.0
Wv = sigmaO2 * PvO2
phi = sigmaCO2*RT*(1.0 + Kc)
phi1 = sigmaCO2*(1+Kc)

PaO2 = np.linspace(40, 150, 100)
VQO2 = (RT / (PiO2 - PaO2)) * (sigmaO2 * (PaO2 - PvO2) + 4 * Z0 * (hemo(sigmaO2 * PaO2) - hemo(Wv)))
PaCO2 = phi * PvCO2 / (VQO2 + phi)

# Plot the RER curve
plt.plot(PaO2, PaCO2, linewidth=2, label='P_aCO2')

# Plot the separate blood and gas RER curves
Rgas = 0.8
Pag1 = Rgas * (PiO2 - PaO2)
plt.plot(PaO2, Pag1, 'r', linewidth=1)
Rgas = 1.8
Pag2 = Rgas * (PiO2 - PaO2)
plt.plot(PaO2, Pag2, '--r', linewidth=1)

Rblood = 0.8
Pab1 = (sigmaCO2 * (1 + Kc)*PvCO2 - Rblood*(sigmaO2*(PaO2 - PvO2) + 4*Z0*(hemo(sigmaO2*PaO2) - hemo(Wv)))) / (sigmaCO2 * (1 + Kc))
plt.plot(PaO2, Pab1, 'g', linewidth=1)
Rblood = 1.8
Pab2 = (phi1*PvCO2 - Rblood*(sigmaO2*(PaO2 - PvO2) + 4*Z0*(hemo(sigmaO2*PaO2) - hemo(Wv)))) / phi1
plt.plot(PaO2, Pab2, '--g', linewidth=1)

plt.ylim(0, 50)
plt.xlabel(r'$P_{a,O2}$ (mm Hg)')
plt.ylabel(r'$P_{a,CO2}$ (mm Hg)')
plt.show()



