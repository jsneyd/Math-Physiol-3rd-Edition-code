
#   -------------------------------------------------------------------
# 
#   This uses the Python module jitcdde to solve delay differential equations
#   for the extended model of cyclic neutropenia. 
#
#   jitcdde doesn't come with standard Python installations 
#   so you may need to install it. ddeint doesn't seem to work, but I don't know why.
# 
#    For Chapter 13, Section 13.1.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from jitcdde import jitcdde, y, t

gam = 0.04
n = 1
tauR= 2.8
K1 = 0.095
b0 = 8
p0 = 0.8
K2 = 0.36
m = 2
tauN = 3.5
A = 10.2
alpha = 2.4

betaRlagR = b0*(K1**n)/(K1**n + y(0,t-tauR)**n)
betaR = b0*(K1**n)/(K1**n + y(0)**n)
phiNlagN = p0*(K2**m)/(K2**m + y(1,t-tauN)**m)
phiN = p0*(K2**m)/(K2**m + y(1)**m)

dRdt = 2*np.exp(-gam*tauR)*betaRlagR*y(0,t-tauR) - (phiN+betaR)*y(0);
dNdt = A*phiNlagN*y(0,t-tauN) - alpha*y(1);

neutropenia_f = [dRdt, dNdt]

dde = jitcdde(neutropenia_f)

ts = np.linspace(0, 200, 10000)
dde.constant_past( [1.0,0.0], time=0.0 )
dde.step_on_discontinuities()
ys = []
for t in ts:
	ys.append(dde.integrate(t))
    
plt.plot(ts, ys, color='red', linewidth=1)