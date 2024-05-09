#  -----------------------------------------------
#  Code to solve a simple version of the Longtin-Milton model of the pupil
#  light reflex. This does not run under Octave, which has not yet implemented
#  the delay differential equation solver, dde23.
#  
#  For Chapter 19, Section 19.7 of 
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#  Written by James Keener and James Sneyd
#  -----------------------------------------------

import numpy as np
from ddeint import ddeint
import matplotlib.pyplot as plt

def model(Y,t,d):
    x = Y(t)
    xd = Y(t-d)
    capF = lambda s: np.heaviside(s,1)*s
    return gamma*capF(np.log(I*areafun(xd))/phibar) - x

def values_before_zero(t):
    return np.array([10,10])

def areafun(s):
    return lam*(theta**n)/(s**n + theta**n)

lam = 30
n = 7
gamma = 5
I = 10
phibar = 1
theta = 10
tau = 1.2

tt = np.linspace(0, 20, 1000)
yy = ddeint(model, values_before_zero, tt, fargs=(tau,))
plt.plot(tt,areafun(yy[:, 0]), lw=2, label="delay = %.01f" % tau)
plt.plot(tt,yy[:, 0], lw=2, label="delay = %.01f" % tau)
plt.legend()