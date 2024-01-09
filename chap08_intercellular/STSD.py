#   -------------------------------------------------------------------
# 
#    Code to calculate the response of a synapse to repeated stimulation.
# 
#    For Chapter 8, Section 8.1.3 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the function for the right-hand side of the ODE
def rhs(n, t):
    return (1 - n) / tau

# First plot a typical time series
p = 0.5
tau = 1
Ninit = 1
N = [Ninit]
T = [0]
delt = 1

for count in range(1, 8):
    tt = np.linspace((count - 1) * delt, count * delt, 50)
    Ni = odeint(rhs, Ninit, tt)
    Ninit = Ni[-1] * (1 - p)
    N.extend(Ni.flatten())
    T.extend(tt)

plt.figure(1)
plt.plot(T, N, linewidth=2)
plt.xlabel('t/tau')
plt.ylabel('n')
plt.show()

# Now plot the frequency versus mean response
fcount = np.arange(-2, 2.1, 0.1)
freq = 10**fcount
delt = 1/freq
nbar = (delt - p * tau + (p ** 2 * tau) / (p + np.exp(delt / tau) - 1)) / delt  # average value
testfit = 1 / (1 + freq * tau * p)

plt.figure(2)
plt.semilogx(freq, nbar, freq, testfit, linewidth=2)
plt.ylim([0, 1])
plt.xlabel('f*tau')
plt.show()
