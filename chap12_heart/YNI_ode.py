
#   -------------------------------------------------------------------
# 
#    ODE integrator for the YNI equations.
# 
#    For Chapter 12, Section 12.2.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define the YNI model
def YNI(t, S):
    V, f, h, p, q, d = S

    # Calculate gating variables
    ah = 0.001209 * np.exp(-(V + 20.) / 6.534)
    bh = 1. / (np.exp(-0.1 * (V + 30.)) + 1.)
    ap = 0.009 / (1. + np.exp(-(V + 3.8) / 9.71)) + 0.0006
    bp = 0.000225 * (V + 40.) / (np.exp((V + 40.) / 13.3) - 1.)
    aq = 0.00034 * (V + 100.) / (np.exp((V + 100.) / 4.4) - 1.) + .0000495
    bq = 0.0005 * (V + 40.) / (1 - np.exp(-V / 6 - 6.6667)) + .0000845
    tm = V + 35.
    ad = (tm > 0) * tm / (1. - np.exp(-tm / 2.5)) + (tm < 0) * tm * np.exp(tm / 2.5) / (
                np.exp(tm / 2.5) - 1.) + (tm == 0) * 2.5
    ad = (V>0)*(.01045*ad + .03125*V/(1-np.exp(-V/4.8))) + \
        (V<0)* 0.01045*ad + .03125*V*np.exp(V/4.8)/(np.exp(V/4.8)-1.) + \
            (V==0)*( .01045*ad +.03125*4.8)
 
    bd = .00421 * (V - 5.) / (np.exp(V / 2.5 - 2.) - 1.)
    af = .000355 * (V + 20.) / (np.exp((V + 20.) / 5.633) - 1.)
    bf = .000944 * (V + 60.) / (1. + np.exp(-(V + 29.5) / 4.16))

    # Calculate gating variable rates
    Fs = af * (1. - f) - bf * f
    Fh = ah * (1. - h) - bh * h
    Fp = ap * (1. - p) - bp * p
    Fq = aq * (1. - q) - bq * q
    Fd = ad * (1. - d) - bd * d

    # Calculate currents
    currents = IV(S)
    INa = currents[0]
    IK = currents[1]
    Ih = currents[2]
    Is = currents[3]
    Il = currents[4]

    # Calculate membrane potential rate
    Fv = -(Is + INa + IK + Ih + Il)

    return [Fv, Fs, Fh, Fp, Fq, Fd]

# Define the ionic currents
def IV(S):
    V, f, h, p, q, d = S

    tm = V + 37.
    am = (tm > 0) * tm / (1. - np.exp(-tm / 10.)) + (tm < 0) * tm * np.exp(tm / 10.) / (
                np.exp(tm / 10.) - 1.) + (tm == 0) * 10
    bm = 40. * np.exp(-0.056 * (V + 62.))
    minf = am / (am + bm)
    INa = gnabar * minf ** 3 * h * (V - Vna)
    Ik = 0.7 * p * (np.exp(0.0277 * (V + 90.)) - 1.) / (np.exp(0.0277 * (V + 40.)))
    Ih = 0.4 * q * (V + 25.)
    Is = 12.5 * (0.95 * f + 0.05) * (np.exp(V / 15. - 2.) - 1.) * (0.95 * d + 0.05)
    Il = 0.8 * (1. - np.exp(-V / 20. - 3.))

    return [INa, Ik, Ih, Is, Il]

# Global parameters
gnabar = 0.5
Vna = 30.

# Initial conditions
V = -43.2572
f = 1
h = 0.9207
p = 0.1003
q = 0.0365
d = 0.0001

# Time parameters
tstep = 1
t_end = 1000
tspan = [0, t_end]

s0 = [V,f,h,p,q,d]

# Solve the ODE
sol = solve_ivp(YNI, tspan, s0, t_eval=np.arange(0, t_end + tstep, tstep))

# Plot results
plt.figure(figsize=(10, 10))
plt.subplot(2, 1, 1)
plt.plot(sol.t, sol.y[0], 'b', label='V')
plt.legend()
plt.xlabel('Time (ms)')
plt.ylabel('V (mV)')
plt.subplot(2, 1, 2)
plt.plot(sol.t, sol.y[1], 'r', label='f')
plt.plot(sol.t, sol.y[2], 'g', label='h')
plt.plot(sol.t, sol.y[3], 'b', label='p')
plt.plot(sol.t, sol.y[4], 'm', label='q')
plt.plot(sol.t, sol.y[5], 'k', label='d')
plt.legend()
plt.xlabel('Time (ms)')
plt.ylabel('Gating Variables')

