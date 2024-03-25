
#   -------------------------------------------------------------------
# 
#   The compartmental model for oxygen uptake by the lungs.
# 
#    For Chapter 14, Section 14.4.5 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fsolve

def Tode(TO, t):
    O, Otilde = get_O_Otilde(TO)
    return D * (sigma * PO - Otilde) - gamma * O


def get_O_Otilde(TO):
    O = fsolve(lambda x: conc(x) - TO, 20)
    Otilde = fsolve(lambda x: D * (sigma * PO - x) + TO - conc(x), 20)
    return O[0], Otilde[0]

def conc(O):
    return O + 4 * THb * O**4 / (KO**4 + O**4)


# Define global parameters
sigma = 1.4
KO = 30 * sigma
THb = 2200
RT = 1.7e-2  # units of uM
D = 32.5
PO = 150  # air partial pressure
gamma = 37

# Solve the ODE for TO
init = 6000
tspan = np.linspace(0, 10, 100)


sol = odeint(Tode, init, tspan)
O = []
Otilde = []
for TO in sol:
    O_val, Otilde_val = get_O_Otilde(TO)
    O.append(O_val)
    Otilde.append(Otilde_val)

# Plot the results
plt.figure(1)
plt.plot(tspan, sol)
plt.xlabel('Time (nondimensional)')
plt.ylabel('[O_2] in Body')

# Calculate and plot oxygen partial pressures
Oa = 104 * sigma
O_values = np.array(O) / sigma
Otilde_values = np.array(Otilde) / sigma

plt.figure(2)
plt.plot(tspan, O_values, 'r', label='O')
plt.plot(tspan, Otilde_values, 'b', label='Otilde')
plt.title('O_2 Partial Pressures')
plt.xlabel('Time (dimensionless)')
plt.legend()
plt.show()

print(f'O_final: {O[-1] / sigma}')
print(f'Otilde_final: {Otilde[-1] / sigma}')
