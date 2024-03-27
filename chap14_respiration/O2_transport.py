
#  Keener and Sneyd. Mathematical Physiology, 3rd edition, Chapter 14,
#  Figure 1.4.

# Written by James Keener and James Sneyd.

#  The graph should be computed in concentrations and then converted to
#  pressures using the solubility. Be careful with the units. It is not a good idea
#  to work in units of Molar. The values for sigma are much too
#  small, and the numerics start going a little weird. So   all
#  units are converted to microMolar. The code plays much more nicely.

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def rhs(x, O):
    fprime = 4 * O[0] ** 3 * KO ** 4 / ((KO ** 4 + O[0] ** 4) ** 2)
    return (PO2 - O) / (1 + 4 * THb * fprime)


def twodrhs(x, dum):
    O, Y = dum
    X = THb - Y
    return [PO2 - O + 4*km1*Y - 4*k1*X*O**4, k1*X*O**4 - km1*Y]


# The figure in the text
sigmaO = 1.4  # Units of uMolar/mm Hg
KO = 30 * sigmaO
PO2 = 104 * sigmaO  # convert pressure to conc (units of uMolar)
THb = 2000  # 2 mM

# The figure in the text
L = 100
xspan = np.linspace(0, L, 2000)
init = [40 * sigmaO]
sol = solve_ivp(rhs, [0, L], init, t_eval=xspan)

plt.figure(1)
plt.plot(sol.t, sol.y[0] / sigmaO, 'r')  # plot the pressure, not the conc
plt.xlabel('dimensionless distance')
plt.ylabel(r'$[{\rm O}_2]$')


# This next part is for one of the Exercises.
# Find the solution without the equilibrium approximation. When solving
# this way, be careful with  initial conditions, and be sure 
# to start on the critical manifold. Otherwise, the solution takes some
# time to catch up and there is an offset in the solution. 

timescale = 10
k1 = timescale
km1 = (KO**4)*k1
Oinit = 40 * sigmaO
Yinit = THb * Oinit ** 4 / (KO**4 + Oinit**4)
init = [Oinit, Yinit]

L = 100
xspan = np.linspace(0, L, 200)
sol = solve_ivp(twodrhs, [0, L], init, method='Radau',t_eval=xspan)

plt.figure(1)
plt.plot(sol.t, sol.y[0] / sigmaO, '--b')  # plot the pressure, not the conc
plt.legend(['With equilibrium approximation', 'Without equilibrium approximation'])
plt.show()



