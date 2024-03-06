
#   -------------------------------------------------------------------
# 
#   This uses the Python module jitcdde to solve delay differential equations
#   for cyclic neutropenia. jitcdde doesn't come with standard Python installations 
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

gam = 0.33
n = 4
tau = 1.39
del_ = 0.08
K1 = 3.07
b0 = 1.62

f1 = b0 * K1 ** n / (K1 ** n + y(0,t-tau) ** n)
f2 = b0 * K1 ** n / (K1 ** n + y(0) ** n)
f = [2 * f1 * y(0,t-tau) * np.exp(-gam * tau) - (del_ + f2) * y(0)]


dde = jitcdde(f)

ts = np.linspace(0, 75, 10000)
dde.constant_past([4.7])
ys = []
for t in ts:
	ys.append(dde.integrate(t))
    
plt.plot(ts, ys, color='red', linewidth=1)