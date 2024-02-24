
#   -------------------------------------------------------------------
# 
#    A simple integrative model, combining multiple CRUs to get a whole-cell
#    response.
# 
#    For Chapter 12, Section 12.3.3 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def Ltrhs(x, t):
    op, n3, cd, Obar, c = x
    alpha = 2 * np.exp(0.012 * (Vm - 35))
    beta = 0.0882 * np.exp(-0.05 * (Vm - 35))
    gamma = 0.44 * cd
    fbar = alpha * f / (alpha + beta)
    omegabar = omega * (beta / 2 + alpha) / (2 * alpha + beta / 2)
    gammabar = gamma * (beta + 2 * alpha) / (alpha + beta)

    # Ca fluxes
    JLtype = Ltype(op)
    Jpm = Vc * c ** 2 / (Kc ** 2 + c ** 2)
    JRyR = RyR_flux(Obar, cd)

    # L-type channel
    dopdt = -g * op + fbar * n3
    dn3dt = -(fbar + gammabar) * n3 + omegabar * (1 - n3 - op) + g * op

    # cd equation
    dcddt = -JLtype + JRyR - D * (cd - c)

    # RyR equation
    dObardt = km2 * (1 - Obar) - k2 * cd * Obar

    # cytosolic calcium equation
    dcdt = Vf * D * (cd - c) - Jpm

    return [dopdt, dn3dt, dcddt, dObardt, dcdt]

def Ltype(op):
    return kL * op * (Vm - 60)

def RyR_flux(Obar, cd):
    return kRyR * Obar * (cd ** 2) / (K1 ** 2 + cd ** 2)



VGCC_mean = 0.05
VGCC_sig = 0.04
kRyR_mean = 20
kRyR_sig = 15

K1 = 0.55
k2 = 0.05
km2 = 0.01
Vc = 0.2
Kc = 0.2
D = 10
Vf = 0.01
f = 0.85
omega = 0.02
g = 2


# Select VGCC and kRyR from a
# uniform distribution and then calculate the responses for all Vstim.
# Do this Nsample times.

Nsample = 100
n = 150
Vstim = np.linspace(-74, 65, n)
peak_JRyR = np.zeros((Nsample, n))
peak_Ltype = np.zeros((Nsample, n))

for j in range(Nsample):
    # first find the unstimulated steady state (at Vm = -75)
    Vm = -75
    kL = 2 * VGCC_sig * np.random.rand() + VGCC_mean - VGCC_sig
    kRyR = 2 * kRyR_sig * np.random.rand() + kRyR_mean - kRyR_sig
    
    dt = 0.01
    tend = 250
    tspan = np.arange(0, tend + dt, dt)
    init = [0.0041, 0.4, 0.006469, 0.75559, 0.006]
    sol = odeint(Ltrhs, init, tspan)
    init = sol[-1, :]      # the steady state for the all the following Vstims

    # Now start at the steady state and stimulate by increasing Vm
    for i in range(n):
        Vm = Vstim[i]
        dt = 0.01
        tend = 10
        tspan = np.arange(0, tend + dt, dt)
        sol = odeint(Ltrhs, init, tspan)
        op, n3, cd, Obar, c = sol.T

        # pull out the peak fluxes for plotting
        peak_JRyR[j, i] = RyR_flux(Obar, cd).max()
        dum = -Ltype(op)
        peak_Ltype[j, i] = dum.max()
        
peak_JRyR_mean = peak_JRyR.mean(axis=0)
peak_Ltype_mean = peak_Ltype.mean(axis=0)    
 
plt.figure(1)
plt.plot(Vstim, peak_JRyR_mean, label='Peak RyR Flux')
plt.plot(Vstim, peak_Ltype_mean, label='Peak L-Type Flux')

 
# Finally, plot the Vstim curve for the parameters in the middle of the
# distribution. Easiest, but less efficient, to recalculate it.


n = 150
peak_JRyR_fix = np.zeros((n))
peak_Ltype_fix = np.zeros((n))
Vstim = np.linspace(-74, 65, n)

# first find the unstimulated steady state (at Vm = -75)
Vm = -75
kL = 0.05
kRyR =  20

dt = 0.01
tend = 250
tspan = np.arange(0, tend + dt, dt)
init = [0.0041, 0.4, 0.006469, 0.75559, 0.006]
sol = odeint(Ltrhs, init, tspan)
init = sol[-1, :]      # the steady state for the all the following Vstims

# Now start at the steady state and stimulate by increasing Vm
for i in range(n):
    Vm = Vstim[i]
    dt = 0.01
    tend = 10
    tspan = np.arange(0, tend + dt, dt)
    sol = odeint(Ltrhs, init, tspan)
    op, n3, cd, Obar, c = sol.T

    # pull out the peak fluxes for plotting
    peak_JRyR_fix[i] = RyR_flux(Obar, cd).max()
    dum = -Ltype(op)
    peak_Ltype_fix[i] = dum.max()
       
plt.plot(Vstim, peak_JRyR_fix, label='Peak RyR (single CRU)')
plt.plot(Vstim, peak_Ltype_fix, label='Peak L-Type (single CRU)')
plt.legend()
plt.xlabel('Stimulating Voltage Clamp Levels')
plt.ylabel('Peak Fluxes')
plt.show()

